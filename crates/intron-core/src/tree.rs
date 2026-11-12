//! Newick trees: parsing, tip lookup, MRCA and re-rooting.
//!
//! Dollo and maximum-parsimony scoring both need a rooted topology, and the root is a *hypothesis*
//! that changes the answer — the same site can need one gain under one rooting and two under
//! another. So [`Tree`] keeps the parsed topology as given, exposes [`Tree::looks_unrooted`], and
//! makes the rooting choice explicit through [`Tree::with_outgroups`] rather than burying it.

use std::io::BufRead;
use std::path::Path;

use rustc_hash::FxHashMap;

use crate::error::tree::TreeError;
use crate::error::{Error, Result};

/// Index of a node in [`Tree::nodes`]. The root is index 0.
pub type NodeId = usize;

#[derive(Clone, Debug)]
pub struct TreeNode {
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
    /// Tip or clade label; bootstrap values land here too.
    pub label: Option<String>,
    /// Length of the edge from this node to its parent.
    pub branch_length: Option<f64>,
}

impl TreeNode {
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

/// A rooted tree held in a flat arena, so traversal is index arithmetic rather than pointers.
#[derive(Clone, Debug)]
pub struct Tree {
    nodes: Vec<TreeNode>,
    tip_index: FxHashMap<String, NodeId>,
    postorder: Vec<NodeId>,
    /// Edges from the root to each node, so `mrca` never recomputes them.
    depth: Vec<u32>,
}

impl Tree {
    pub fn root(&self) -> NodeId {
        0
    }

    pub fn node(&self, id: NodeId) -> &TreeNode {
        &self.nodes[id]
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn children(&self, id: NodeId) -> &[NodeId] {
        &self.nodes[id].children
    }

    pub fn parent(&self, id: NodeId) -> Option<NodeId> {
        self.nodes[id].parent
    }

    pub fn label(&self, id: NodeId) -> Option<&str> {
        self.nodes[id].label.as_deref()
    }

    pub fn tips(&self) -> Vec<NodeId> {
        (0..self.nodes.len())
            .filter(|&i| self.nodes[i].is_leaf())
            .collect()
    }

    pub fn n_tips(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf()).count()
    }

    /// Tip ids in left-to-right (Newick) order.
    ///
    /// A postorder lists children before parents and siblings left-to-right, so filtering it
    /// forward yields the tips in display order.
    pub fn tip_labels(&self) -> Vec<String> {
        self.postorder
            .iter()
            .copied()
            .filter(|&i| self.nodes[i].is_leaf())
            .map(|i| self.nodes[i].label.clone().unwrap_or_default())
            .collect()
    }

    pub fn tip_index(&self, label: &str) -> Option<NodeId> {
        self.tip_index.get(label).copied()
    }

    /// Nodes in postorder (children before parents). The root is last.
    pub fn postorder(&self) -> &[NodeId] {
        &self.postorder
    }

    /// Every edge, identified by its child node; the root has no incoming edge.
    pub fn edges(&self) -> impl Iterator<Item = (NodeId, NodeId)> + '_ {
        (1..self.nodes.len()).filter_map(move |c| self.nodes[c].parent.map(|p| (p, c)))
    }

    /// `true` when the root has three or more children, i.e. the Newick encodes an unrooted tree.
    pub fn looks_unrooted(&self) -> bool {
        self.nodes[0].children.len() >= 3
    }

    /// All leaf descendants of `node`.
    pub fn leaves_under(&self, node: NodeId) -> Vec<NodeId> {
        let mut out = Vec::new();
        let mut stack = vec![node];
        while let Some(n) = stack.pop() {
            if self.nodes[n].is_leaf() {
                out.push(n);
            } else {
                stack.extend(self.nodes[n].children.iter().rev().copied());
            }
        }
        out.sort_unstable();
        out
    }

    /// Most recent common ancestor of a set of nodes: the set is folded through a two-node ancestor
    /// search over the precomputed depths, O(k · height).
    pub fn mrca(&self, nodes: &[NodeId]) -> Option<NodeId> {
        let mut it = nodes.iter().copied();
        let mut acc = it.next()?;
        for n in it {
            acc = self.lca(acc, n)?;
        }
        Some(acc)
    }

    fn lca(&self, mut a: NodeId, mut b: NodeId) -> Option<NodeId> {
        while self.depth[a] > self.depth[b] {
            a = self.nodes[a].parent?;
        }
        while self.depth[b] > self.depth[a] {
            b = self.nodes[b].parent?;
        }
        while a != b {
            a = self.nodes[a].parent?;
            b = self.nodes[b].parent?;
        }
        Some(a)
    }

    /// One readable, unique name per node, for tables.
    ///
    /// The root is `root`; a tip is its label; an internal node is its non-numeric label (numeric ones
    /// are support values) or else its tip list. Names that would collide — two clades printing the same
    /// truncated tip list, or an internal node that is itself labelled `root` — get `#<node id>` appended,
    /// so a name always identifies exactly one node.
    pub fn node_names(&self) -> Vec<String> {
        let mut names: Vec<String> = (0..self.len())
            .map(|v| {
                if v == self.root() {
                    "root".to_string()
                } else {
                    self.edge_name(v)
                }
            })
            .collect();
        let mut count: FxHashMap<String, u32> = FxHashMap::default();
        for n in &names {
            *count.entry(n.clone()).or_default() += 1;
        }
        for (v, n) in names.iter_mut().enumerate() {
            if v != self.root() && (count[n.as_str()] > 1 || n == "root") {
                n.push_str(&format!("#{v}"));
            }
        }
        names
    }

    /// `Some(mrca)` when exactly `nodes` form a clade of the tree, `None` otherwise.
    ///
    /// This is the quick Dollo pre-test: a site gained once and never regained must define a clade.
    pub fn monophyletic_group(&self, nodes: &[NodeId]) -> Option<NodeId> {
        let m = self.mrca(nodes)?;
        let under = self.leaves_under(m);
        (under.len() == nodes.len()).then_some(m)
    }

    /// Human-readable name for an edge: the child's label, or its descendant tip list.
    pub fn edge_name(&self, child: NodeId) -> String {
        // A numeric label on an internal node is a bootstrap/support value, not a name: printing
        // `100` as the branch where an intron was gained would be ambiguous and meaningless.
        if let Some(l) = self
            .label(child)
            .filter(|l| self.node(child).is_leaf() || l.parse::<f64>().is_err())
        {
            return l.to_string();
        }
        let tips = self.leaves_under(child);
        let names: Vec<String> = tips
            .iter()
            .take(4)
            .map(|&t| self.label(t).unwrap_or("?").to_string())
            .collect();
        if tips.len() > 4 {
            format!("{}…({} tips)", names.join(","), tips.len())
        } else {
            names.join(",")
        }
    }

    /// Re-root on the undirected edge between two adjacent nodes.
    ///
    /// The edge is split by a brand-new root, so every original tip stays a tip — rooting *at* an
    /// outgroup tip would otherwise turn that tip into an internal node and silently drop it from the
    /// analysis. The two halves inherit half the original branch length each, which preserves
    /// pairwise path lengths.
    pub fn reroot_on_edge(&self, a: NodeId, b: NodeId) -> Result<Tree> {
        if !self.neighbours(a).contains(&b) {
            return Err(Error::Tree(TreeError::Unrootable(format!(
                "nodes {a} and {b} are not joined by an edge"
            ))));
        }
        let half = self.edge_length(a, b).map(|l| l / 2.0);
        let mut nodes = vec![TreeNode {
            parent: None,
            children: Vec::new(),
            label: None,
            branch_length: None,
        }];
        let ai = self.grow(a, Some(b), 0, half, &mut nodes);
        let bi = self.grow(b, Some(a), 0, half, &mut nodes);
        nodes[0].children = vec![ai, bi];
        let mut t = Tree {
            nodes,
            tip_index: Default::default(),
            postorder: Vec::new(),
            depth: Vec::new(),
        };
        t.rebuild();
        Ok(t)
    }

    /// Re-root on the branch immediately above `node`.
    pub fn reroot_above(&self, node: NodeId) -> Result<Tree> {
        let p = self.parent(node).ok_or_else(|| {
            Error::Tree(TreeError::Unrootable(
                "the requested node is already the root".into(),
            ))
        })?;
        self.reroot_on_edge(node, p)
    }

    fn neighbours(&self, n: NodeId) -> Vec<NodeId> {
        let mut v = self.nodes[n].children.clone();
        if let Some(p) = self.nodes[n].parent {
            v.push(p);
        }
        v
    }

    /// Length of the edge joining two adjacent nodes, stored on whichever end is the child.
    fn edge_length(&self, a: NodeId, b: NodeId) -> Option<f64> {
        if self.nodes[b].parent == Some(a) {
            self.nodes[b].branch_length
        } else if self.nodes[a].parent == Some(b) {
            self.nodes[a].branch_length
        } else {
            None
        }
    }

    /// Copy the component of the tree reached from `src` without crossing `exclude`.
    fn grow(
        &self,
        src: NodeId,
        exclude: Option<NodeId>,
        parent: NodeId,
        edge_len: Option<f64>,
        out: &mut Vec<TreeNode>,
    ) -> NodeId {
        let onward: Vec<NodeId> = self
            .neighbours(src)
            .into_iter()
            .filter(|n| Some(*n) != exclude)
            .collect();
        // Re-rooting a rooted binary tree leaves its old root with a single continuing neighbour. Keep
        // it and the same tip set would appear as two nested "clades", so splice it out and merge the
        // two branch lengths.
        if onward.len() == 1 {
            let nb = onward[0];
            let l = if self.nodes[nb].parent == Some(src) {
                self.nodes[nb].branch_length
            } else {
                self.nodes[src].branch_length
            };
            let merged = edge_len.zip(l).map(|(x, y)| x + y).or(edge_len).or(l);
            return self.grow(nb, Some(src), parent, merged, out);
        }
        let id = out.len();
        out.push(TreeNode {
            parent: Some(parent),
            children: Vec::new(),
            label: self.nodes[src].label.clone(),
            branch_length: edge_len,
        });
        for nb in onward {
            // Travelling down uses the child's length; travelling up crosses src's own edge.
            let l = if self.nodes[nb].parent == Some(src) {
                self.nodes[nb].branch_length
            } else {
                self.nodes[src].branch_length
            };
            let c = self.grow(nb, Some(src), id, l, out);
            out[id].children.push(c);
        }
        id
    }

    /// Re-root so the given outgroup labels are outside everything else.
    ///
    /// Accepts one tip (the common case) or a monophyletic set, in which case the branch above their
    /// MRCA is used. Errors if any label is not a tip or the labels do not form one clade.
    pub fn with_outgroups(&self, outgroups: &[String]) -> Result<Tree> {
        if outgroups.is_empty() {
            return Ok(self.clone());
        }
        let mut nodes = Vec::with_capacity(outgroups.len());
        for name in outgroups {
            let Some(n) = self.tip_index(name) else {
                return Err(Error::Tree(TreeError::UnknownOutgroup(name.clone())));
            };
            nodes.push(n);
        }
        let target = if nodes.len() == 1 {
            nodes[0]
        } else {
            let m = self.mrca(&nodes).ok_or_else(|| {
                Error::Tree(TreeError::Unrootable(
                    "outgroup taxa have no common ancestor".into(),
                ))
            })?;
            let under = self.leaves_under(m);
            if under.len() == nodes.len() {
                m
            } else {
                // On an unrooted tree (trifurcating root) the outgroup may be a valid split without
                // being a clade: {A,B} | {C,D} in `(A,B,(C,D));`. Then its complement is the clade.
                let outs: std::collections::HashSet<NodeId> = nodes.iter().copied().collect();
                let rest: Vec<NodeId> = self
                    .tips()
                    .into_iter()
                    .filter(|t| !outs.contains(t))
                    .collect();
                let comp = if rest.is_empty() {
                    None
                } else {
                    self.mrca(&rest)
                };
                match comp {
                    Some(c)
                        if self.leaves_under(c).len() == rest.len() && self.parent(c).is_some() =>
                    {
                        c
                    }
                    _ => {
                        return Err(Error::Tree(TreeError::Unrootable(format!(
                            "outgroup is neither a clade nor the complement of one: its MRCA also contains {} other tip(s)",
                            under.len() - nodes.len()
                        ))));
                    }
                }
            }
        };
        self.reroot_above(target)
    }

    fn rebuild(&mut self) {
        self.postorder = default_postorder(&self.nodes);
        // Parents precede children in the reversed postorder, so one pass fills every depth.
        let mut depth = vec![0u32; self.nodes.len()];
        for &v in self.postorder.iter().rev() {
            if let Some(p) = self.nodes[v].parent {
                depth[v] = depth[p] + 1;
            }
        }
        self.depth = depth;
        self.tip_index = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.is_leaf())
            .filter_map(|(i, n)| n.label.clone().map(|l| (l, i)))
            .collect();
    }

    /// Renumber so the root is index 0, keeping labels and branch lengths attached to nodes.
    fn move_root_to_front(&mut self, old_root: NodeId) {
        // Build a fresh arena by preorder traversal from `old_root`; the root is then index 0.
        let mut map = vec![usize::MAX; self.nodes.len()];
        let mut nodes: Vec<TreeNode> = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![(old_root, None::<usize>)];
        while let Some((src, parent)) = stack.pop() {
            let dst = nodes.len();
            map[src] = dst;
            nodes.push(TreeNode {
                parent,
                children: Vec::new(),
                label: self.nodes[src].label.clone(),
                branch_length: self.nodes[src].branch_length,
            });
            if let Some(p) = parent {
                nodes[p].children.push(dst);
            }
            for &c in self.nodes[src].children.iter().rev() {
                stack.push((c, Some(dst)));
            }
        }
        self.nodes = nodes;
        self.rebuild();
    }

    /// Serialise back to Newick with branch lengths and internal labels.
    pub fn to_newick(&self) -> String {
        let mut out = String::new();
        write_node(self, self.root(), &mut out);
        out.push(';');
        out
    }

    /// Read a Newick file, ignoring blank lines, comments and any trailing trees.
    pub fn read(path: impl AsRef<Path>) -> Result<Tree> {
        let path = path.as_ref();
        let mut reader = crate::util::open_buffered(path, 1 << 16)?;
        let mut text = String::new();
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let t = line.trim();
            // A line that is only a comment (`[ ... ]`) is skipped; a `[&R] (...);` tree line is not.
            let comment_only = t.starts_with('[') && t.ends_with(']') && !t.contains('(');
            if t.is_empty() || t.starts_with('#') || comment_only {
                continue;
            }
            text.push_str(t);
            if t.ends_with(';') {
                break;
            }
        }
        let tree = Tree::from_newick(&text).map_err(|e| match e {
            TreeError::Parse {
                offset,
                message,
                snippet,
            } => TreeError::Parse {
                offset,
                message,
                snippet: format!("{snippet} (in {})", path.display()),
            },
            other => other,
        })?;
        tree.validate()?;
        Ok(tree)
    }

    /// Parse a Newick string.
    ///
    /// ```
    /// use intron_core::tree::Tree;
    /// let t = Tree::from_newick("((human,chimp)homo,(gorilla)gorilla)root;").unwrap();
    /// assert_eq!(t.n_tips(), 3);
    /// assert_eq!(t.tip_index("chimp"), Some(t.tip_index("chimp").unwrap()));
    /// let m = t.mrca(&[t.tip_index("human").unwrap(), t.tip_index("chimp").unwrap()]).unwrap();
    /// assert_eq!(t.label(m), Some("homo"));
    /// assert_eq!(t.monophyletic_group(&[t.tip_index("human").unwrap(), t.tip_index("chimp").unwrap()]), Some(m));
    /// ```
    pub fn from_newick(s: &str) -> std::result::Result<Tree, TreeError> {
        let bytes = s.as_bytes();
        let mut p = Parser { b: bytes, i: 0 };
        p.skip_ws();
        if p.i >= bytes.len() {
            return Err(TreeError::Parse {
                offset: 0,
                message: "empty tree".into(),
                snippet: s.into(),
            });
        }
        let mut nodes = Vec::new();
        let root = p.parse_node(&mut nodes)?;
        p.skip_ws();
        // A trailing ';' is conventional; anything else after the first tree is ignored.
        match nodes.len() {
            0 => return Err(TreeError::Empty),
            1 => {
                return Err(TreeError::Parse {
                    offset: p.i,
                    message: "a single tip is not a tree".into(),
                    snippet: snippet(s, p.i),
                })
            }
            _ => {}
        }
        let mut t = Tree {
            nodes,
            tip_index: Default::default(),
            postorder: Vec::new(),
            depth: Vec::new(),
        };
        if root != 0 {
            // The parser appends the root last; rotate so the arena convention holds.
            t.move_root_to_front(root);
        }
        t.rebuild();
        t.validate()?;
        Ok(t)
    }

    /// Structural checks that make a bad tree fail loudly instead of scoring silently.
    fn validate(&self) -> std::result::Result<(), TreeError> {
        let tips = self.n_tips();
        if tips == 0 {
            return Err(TreeError::Empty);
        }
        if self.tip_index.len() != tips {
            let mut seen = std::collections::HashSet::new();
            for t in self.tips() {
                if let Some(l) = self.label(t) {
                    if !seen.insert(l.to_string()) {
                        return Err(TreeError::DuplicateTip(l.to_string()));
                    }
                }
            }
            return Err(TreeError::DuplicateTip(String::new()));
        }
        Ok(())
    }
}

fn default_postorder(nodes: &[TreeNode]) -> Vec<NodeId> {
    // Children before parents, discovered by an explicit stack from the root (index 0 or parent=None).
    let root = nodes.iter().position(|n| n.parent.is_none()).unwrap_or(0);
    let mut order = Vec::with_capacity(nodes.len());
    let mut seen = vec![false; nodes.len()];
    let mut stack = vec![(root, false)];
    while let Some((n, expanded)) = stack.pop() {
        if expanded {
            order.push(n);
            continue;
        }
        if seen[n] {
            continue;
        }
        seen[n] = true;
        stack.push((n, true));
        for &c in nodes[n].children.iter().rev() {
            stack.push((c, false));
        }
    }
    order
}

fn write_node(t: &Tree, n: NodeId, out: &mut String) {
    let node = t.node(n);
    if !node.children.is_empty() {
        out.push('(');
        for (i, &c) in node.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            write_node(t, c, out);
        }
        out.push(')');
    }
    if let Some(l) = &node.label {
        if l.bytes()
            .any(|b| matches!(b, b'(' | b')' | b',' | b':' | b';' | b' ' | b'\'' | b'_'))
        {
            out.push('\'');
            out.push_str(&l.replace('\'', "''"));
            out.push('\'');
        } else {
            out.push_str(l);
        }
    }
    if let Some(bl) = node.branch_length {
        let _ = write_bl(out, bl);
    }
}

fn write_bl(out: &mut String, bl: f64) -> std::fmt::Result {
    use std::fmt::Write;
    write!(out, ":{bl}")
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Parser<'a> {
    /// Skip whitespace and Newick/NHX comments (`[...]`, which may nest). Comments carry annotations
    /// such as `[&R]`, `[&rate=1]` or `[&&NHX:S=x]`; none of them is part of a label or a length.
    fn skip_ws(&mut self) {
        loop {
            while self.i < self.b.len() && (self.b[self.i] as char).is_whitespace() {
                self.i += 1;
            }
            if self.peek() == Some(b'[') {
                let mut depth = 0usize;
                while let Some(c) = self.peek() {
                    self.i += 1;
                    match c {
                        b'[' => depth += 1,
                        b']' => {
                            depth -= 1;
                            if depth == 0 {
                                break;
                            }
                        }
                        _ => {}
                    }
                }
            } else {
                break;
            }
        }
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.i).copied()
    }

    fn err(&self, msg: impl Into<String>) -> TreeError {
        TreeError::Parse {
            offset: self.i,
            message: msg.into(),
            snippet: snippet(
                std::str::from_utf8(self.b).unwrap_or("<invalid utf8>"),
                self.i,
            ),
        }
    }

    /// Parse one node (with its descendants), returning its arena index.
    fn parse_node(&mut self, nodes: &mut Vec<TreeNode>) -> std::result::Result<NodeId, TreeError> {
        self.skip_ws();
        let mut children = Vec::new();
        if self.peek() == Some(b'(') {
            self.i += 1;
            loop {
                let c = self.parse_node(nodes)?;
                children.push(c);
                self.skip_ws();
                match self.peek() {
                    Some(b',') => self.i += 1,
                    Some(b')') => {
                        self.i += 1;
                        break;
                    }
                    _ => return Err(self.err("expected ',' or ')' inside a clade")),
                }
            }
        }
        let label = self.parse_label()?;
        self.skip_ws(); // a comment may sit between the label and `:length`
        let branch_length = self.parse_branch_length()?;
        let id = nodes.len();
        nodes.push(TreeNode {
            parent: None,
            children,
            label,
            branch_length,
        });
        // Two separate statements: the read of `children[i]` ends before the write to `nodes[c]`.
        for i in 0..nodes[id].children.len() {
            let c = nodes[id].children[i];
            nodes[c].parent = Some(id);
        }
        Ok(id)
    }

    fn parse_label(&mut self) -> std::result::Result<Option<String>, TreeError> {
        self.skip_ws();
        let start = self.i;
        if self.peek() == Some(b'\'') {
            self.i += 1;
            // Collect bytes and decode once: pushing `byte as char` would turn UTF-8 labels such as
            // `Brésil` into mojibake that no longer matches the matrix's taxon names.
            let mut s: Vec<u8> = Vec::new();
            while let Some(c) = self.peek() {
                self.i += 1;
                if c == b'\'' {
                    if self.peek() == Some(b'\'') {
                        s.push(b'\'');
                        self.i += 1;
                        continue;
                    }
                    return String::from_utf8(s)
                        .map(Some)
                        .map_err(|_| self.err("non-UTF8 quoted label"));
                }
                s.push(c);
            }
            return Err(self.err("unterminated quoted label"));
        }
        while let Some(c) = self.peek() {
            if matches!(c, b'(' | b')' | b',' | b':' | b';' | b'[') {
                break;
            }
            self.i += 1;
        }
        if self.i == start {
            return Ok(None);
        }
        let raw =
            std::str::from_utf8(&self.b[start..self.i]).map_err(|_| self.err("non-UTF8 label"))?;
        Ok(Some(raw.trim().to_string()))
    }

    fn parse_branch_length(&mut self) -> std::result::Result<Option<f64>, TreeError> {
        if self.peek() != Some(b':') {
            return Ok(None);
        }
        self.i += 1;
        let start = self.i;
        while let Some(c) = self.peek() {
            if matches!(c, b'(' | b')' | b',' | b';' | b'[') {
                break;
            }
            self.i += 1;
        }
        let text = std::str::from_utf8(&self.b[start..self.i])
            .map_err(|_| self.err("non-UTF8 branch length"))?;
        let text = text.trim();
        if text.is_empty() {
            return Err(self.err("':' with no branch length"));
        }
        text.parse::<f64>()
            .map(Some)
            .map_err(|_| self.err(format!("bad branch length {text:?}")))
    }
}

fn snippet(s: &str, at: usize) -> String {
    // Move both ends to character boundaries: byte offsets can land inside a multibyte label.
    let mut start = at.saturating_sub(20).min(s.len());
    while !s.is_char_boundary(start) {
        start -= 1;
    }
    let mut end = (at + 20).min(s.len());
    while !s.is_char_boundary(end) {
        end += 1;
    }
    format!("…{}…", &s[start..end])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(t: &Tree) -> Vec<String> {
        t.tip_labels()
    }

    /// The pre-optimisation algorithm, kept as the reference the fast `mrca` must agree with.
    fn naive_mrca(t: &Tree, nodes: &[NodeId]) -> Option<NodeId> {
        fn depth(t: &Tree, mut n: NodeId) -> usize {
            let mut d = 0;
            while let Some(p) = t.parent(n) {
                d += 1;
                n = p;
            }
            d
        }
        let mut cur = nodes.to_vec();
        loop {
            if cur.is_empty() {
                return None;
            }
            if cur.iter().all(|&n| n == cur[0]) {
                return Some(cur[0]);
            }
            let deepest = (0..cur.len()).max_by_key(|&i| depth(t, cur[i])).unwrap();
            cur[deepest] = t.parent(cur[deepest])?;
        }
    }

    #[test]
    fn fast_mrca_agrees_with_the_naive_algorithm_on_random_trees() {
        let mut state = 0x9E3779B97F4A7C15u64;
        let mut next = move |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        for _ in 0..60 {
            // Random tree with polytomies: repeatedly merge a random handful of open clades.
            let n_tips = 3 + next(40);
            let mut open: Vec<String> = (0..n_tips).map(|i| format!("t{i}")).collect();
            while open.len() > 1 {
                let k = (2 + next(3)).min(open.len());
                let mut group = Vec::new();
                for _ in 0..k {
                    let i = next(open.len());
                    group.push(open.swap_remove(i));
                }
                open.push(format!("({})", group.join(",")));
            }
            let nwk = format!("{};", open[0]);
            let Ok(t) = Tree::from_newick(&nwk) else {
                continue;
            };
            let tips = t.tips();
            for _ in 0..40 {
                let k = 1 + next(tips.len().min(12));
                let set: Vec<NodeId> = (0..k).map(|_| tips[next(tips.len())]).collect();
                assert_eq!(t.mrca(&set), naive_mrca(&t, &set), "{nwk} {set:?}");
            }
        }
    }

    #[test]
    fn newick_comments_are_ignored_everywhere() {
        let t = Tree::from_newick("[&R] ((A[&&NHX:S=x]:0.1[&rate=1],B)90,C);").unwrap();
        assert_eq!(ids(&t), vec!["A", "B", "C"]);
        assert_eq!(
            t.edge_length(
                t.tip_index("A").unwrap(),
                t.parent(t.tip_index("A").unwrap()).unwrap()
            ),
            Some(0.1)
        );
    }

    #[test]
    fn quoted_labels_keep_their_utf8() {
        let t = Tree::from_newick("('Brésil',B,C);").unwrap();
        assert!(t.tip_index("Brésil").is_some(), "{:?}", ids(&t));
    }

    #[test]
    fn parse_errors_with_multibyte_text_do_not_panic() {
        let long = "é".repeat(30);
        let err = Tree::from_newick(&format!("(A:1é,Bx{long}")).unwrap_err();
        assert!(err.to_string().contains("newick parse error"), "{err}");
    }

    #[test]
    fn numeric_internal_labels_are_support_not_names() {
        let t = Tree::from_newick("((A,B)100,C);").unwrap();
        let ab = t.parent(t.tip_index("A").unwrap()).unwrap();
        assert_eq!(t.edge_name(ab), "A,B");
        let named = Tree::from_newick("((A,B)ab,C);").unwrap();
        let ab = named.parent(named.tip_index("A").unwrap()).unwrap();
        assert_eq!(named.edge_name(ab), "ab");
    }

    #[test]
    fn multi_taxon_outgroup_may_be_the_complement_of_a_clade() {
        // {A,B} is a valid split of this unrooted tree although it is not a clade.
        let t = Tree::from_newick("(A,B,(C,D));").unwrap();
        let r = t.with_outgroups(&["A".into(), "B".into()]).unwrap();
        let a = r.tip_index("A").unwrap();
        let b = r.tip_index("B").unwrap();
        assert!(r.monophyletic_group(&[a, b]).is_some());
        assert_eq!(r.n_tips(), 4);
        assert!(
            t.with_outgroups(&["A".into(), "C".into()]).is_err(),
            "{{A,C}} is not a split"
        );
    }

    #[test]
    fn rerooting_a_rooted_binary_tree_leaves_no_unary_node() {
        let t = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let r = t.with_outgroups(&["A".into()]).unwrap();
        for v in 0..r.len() {
            assert_ne!(
                r.children(v).len(),
                1,
                "node {v} is unary in {}",
                r.to_newick()
            );
        }
        assert_eq!(r.n_tips(), 4);
    }

    #[test]
    fn parses_balanced_topology_with_labels() {
        let t = Tree::from_newick("((A,B),(C,(D,E))F);").unwrap();
        assert_eq!(t.n_tips(), 5);
        assert_eq!(ids(&t), vec!["A", "B", "C", "D", "E"]);
        assert_eq!(t.len(), 9, "5 tips + 3 internal + root");
        assert_eq!(t.node(t.root()).children.len(), 2);
    }

    #[test]
    fn parses_branch_lengths_and_bootstrap() {
        let t = Tree::from_newick("((A:0.1,B:0.2)85:0.3,C:0.4)root;").unwrap();
        let a = t.tip_index("A").unwrap();
        assert_eq!(t.node(a).branch_length, Some(0.1));
        let clade = t.parent(a).unwrap();
        assert_eq!(t.label(clade), Some("85"));
        assert_eq!(t.node(clade).branch_length, Some(0.3));
        assert_eq!(t.label(t.root()), Some("root"));
    }

    #[test]
    fn handles_quoted_and_odd_labels() {
        // A trailing ';' is conventional, not required.
        let t = Tree::from_newick("(('in group',B_1)x)root").unwrap();
        assert_eq!(t.n_tips(), 2);
        assert!(t.tip_index("in group").is_some());
        assert!(t.tip_index("B_1").is_some());
        assert_eq!(
            t.label(t.parent(t.tip_index("in group").unwrap()).unwrap()),
            Some("x")
        );
    }

    #[test]
    fn tolerates_whitespace_and_newlines() {
        let t = Tree::from_newick(" ( (A , B) , C ) ;\n").unwrap();
        assert_eq!(t.n_tips(), 3);
    }

    #[test]
    fn rejects_malformed_newick() {
        for bad in ["((A,B)", "(A,B", "(A,,B)", ":", "((A:0.1,B:),C);"] {
            assert!(Tree::from_newick(bad).is_err(), "{bad:?} should fail");
        }
        let e = Tree::from_newick("(A,B").unwrap_err();
        assert!(e.to_string().contains("expected"), "{e}");
    }

    #[test]
    fn mrca_and_monophyly() {
        let t = Tree::from_newick("(((A,B),(C,D)),E);").unwrap();
        let (a, b, c, e) = (
            t.tip_index("A").unwrap(),
            t.tip_index("B").unwrap(),
            t.tip_index("C").unwrap(),
            t.tip_index("E").unwrap(),
        );
        let ab = t.mrca(&[a, b]).unwrap();
        assert_eq!(t.leaves_under(ab).len(), 2);
        assert_eq!(t.monophyletic_group(&[a, b]), Some(ab));
        assert_eq!(
            t.monophyletic_group(&[a, c]),
            None,
            "A and C are not a clade by themselves"
        );
        assert_eq!(t.mrca(&[a, b, c]), t.parent(ab));
        assert_eq!(t.monophyletic_group(&[a, b, c]), None);
        assert_eq!(t.mrca(&[]), None);
        assert!(
            t.monophyletic_group(&[e]).is_some(),
            "a single tip is a clade"
        );
    }

    #[test]
    fn unrooted_detection() {
        assert!(Tree::from_newick("(A,B,C);").unwrap().looks_unrooted());
        assert!(!Tree::from_newick("((A,B),C);").unwrap().looks_unrooted());
    }

    #[test]
    fn reroot_above_a_tip_splits_the_branch_and_keeps_every_tip() {
        let t = Tree::from_newick("(((A,B),(C,D)),E);").unwrap();
        let e = t.tip_index("E").unwrap();
        let r = t.reroot_above(e).unwrap();
        assert_eq!(
            r.n_tips(),
            5,
            "the outgroup stays a tip instead of becoming the root"
        );
        let kids = r.children(r.root());
        assert_eq!(kids.len(), 2, "a fresh binary root sits on the split edge");
        assert_eq!(r.leaves_under(kids[0]).len(), 1);
        assert_eq!(r.leaves_under(kids[1]).len(), 4);
        // Ingroup relationships are untouched by where the tree is rooted.
        let (a, b) = (r.tip_index("A").unwrap(), r.tip_index("B").unwrap());
        assert_eq!(r.leaves_under(r.mrca(&[a, b]).unwrap()).len(), 2);
        let (c, d) = (r.tip_index("C").unwrap(), r.tip_index("D").unwrap());
        assert!(
            r.monophyletic_group(&[a, b, c, d]).is_some(),
            "the ingroup is now a clade"
        );
    }

    #[test]
    fn reroot_splits_the_edge_length_evenly_and_leaves_other_edges_alone() {
        let t = Tree::from_newick("(((A:1,B:2)x:3,(C:4,D:5)y:6)z,E:7);").unwrap();
        let e = t.tip_index("E").unwrap();
        let r = t.reroot_above(e).unwrap();
        let kids = r.children(r.root());
        assert_eq!(
            r.node(kids[0]).branch_length,
            Some(3.5),
            "half of E's 7.0 branch"
        );
        assert_eq!(r.node(kids[1]).branch_length, Some(3.5));
        assert_eq!(r.node(r.tip_index("A").unwrap()).branch_length, Some(1.0));
        assert_eq!(r.node(r.tip_index("D").unwrap()).branch_length, Some(5.0));
        assert_eq!(r.tip_labels().len(), 5);
    }

    #[test]
    fn reroot_inside_a_clade_reverses_only_the_path() {
        let t = Tree::from_newick("(((A,B)ab,(C,D)cd)abcd,(E,F)ef)root;").unwrap();
        let c = t.tip_index("C").unwrap();
        let r = t.reroot_above(c).unwrap();
        assert_eq!(r.n_tips(), 6);
        assert_eq!(
            r.tip_labels().iter().filter(|s| s.as_str() == "C").count(),
            1
        );
        let kids = r.children(r.root());
        assert_eq!(
            r.leaves_under(kids[0]).len(),
            1,
            "C alone on one side of the new root"
        );
        assert_eq!(r.leaves_under(kids[1]).len(), 5);
        let (a, b) = (r.tip_index("A").unwrap(), r.tip_index("B").unwrap());
        assert_eq!(
            r.leaves_under(r.mrca(&[a, b]).unwrap()).len(),
            2,
            "A and B are still sisters"
        );
        let (e, f) = (r.tip_index("E").unwrap(), r.tip_index("F").unwrap());
        assert!(r.monophyletic_group(&[e, f]).is_some(), "E,F unaffected");
    }

    #[test]
    fn reroot_on_a_non_edge_is_rejected() {
        let t = Tree::from_newick("(((A,B),(C,D)),E);").unwrap();
        let (a, c) = (t.tip_index("A").unwrap(), t.tip_index("C").unwrap());
        assert!(t.reroot_on_edge(a, c).is_err(), "A and C are not adjacent");
        assert!(
            t.reroot_above(t.root()).is_err(),
            "the root has no branch above it"
        );
    }

    #[test]
    fn with_outgroups_validates_labels_and_monophyly() {
        let t = Tree::from_newick("(((A,B),(C,D)),E);").unwrap();
        let r = t.with_outgroups(&["E".into()]).unwrap();
        let kids = r.children(r.root());
        assert_eq!(
            r.leaves_under(kids[0]),
            vec![r.tip_index("E").unwrap()],
            "E sits alone on one side"
        );
        assert!(t.with_outgroups(&["Z".into()]).is_err(), "unknown outgroup");
        assert!(
            t.with_outgroups(&["A".into(), "C".into()]).is_err(),
            "non-monophyletic outgroup"
        );
        let both = t.with_outgroups(&["A".into(), "B".into()]).unwrap();
        assert_eq!(
            both.leaves_under(both.root()).len(),
            5,
            "the whole tree is preserved"
        );
        assert!(both
            .monophyletic_group(&[both.tip_index("C").unwrap(), both.tip_index("D").unwrap()])
            .is_some());
    }

    #[test]
    fn round_trips_through_newick() {
        let src = "((A:0.10,B:0.20)ab:0.30,(C:0.40,D:0.50)cd:0.60)root;";
        let t = Tree::from_newick(src).unwrap();
        let out = t.to_newick();
        let t2 = Tree::from_newick(&out).unwrap();
        assert_eq!(t2.tip_labels(), t.tip_labels());
        assert_eq!(t2.n_tips(), 4);
        assert_eq!(t2.node(t2.tip_index("A").unwrap()).branch_length, Some(0.1));
    }

    #[test]
    fn edge_names_use_labels_then_tip_sets() {
        let t = Tree::from_newick("(((A,B),(C,D))inner,E);").unwrap();
        let ab = t
            .mrca(&[t.tip_index("A").unwrap(), t.tip_index("B").unwrap()])
            .unwrap();
        assert_eq!(t.edge_name(t.tip_index("A").unwrap()), "A");
        assert_eq!(
            t.edge_name(ab),
            "A,B",
            "unlabelled clades are described by their tips"
        );
        let labelled = t
            .mrca(&[t.tip_index("C").unwrap(), t.tip_index("D").unwrap()])
            .unwrap();
        assert_eq!(t.edge_name(labelled), "C,D");
        let root_name = t.edge_name(t.root());
        assert!(
            root_name.contains("5 tips"),
            "the root spans everything: {root_name}"
        );
    }

    #[test]
    fn reads_a_file_and_reports_the_path_in_errors() {
        let dir = std::env::temp_dir().join(format!("ip_tree_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let ok = dir.join("t.nwk");
        std::fs::write(&ok, "((A,B),C);\n").unwrap();
        let t = Tree::read(&ok).unwrap();
        assert_eq!(t.n_tips(), 3);
        let bad = dir.join("bad.nwk");
        std::fs::write(&bad, "((A,B\n").unwrap();
        let e = Tree::read(&bad).unwrap_err();
        assert!(
            e.to_string().contains("bad.nwk"),
            "error should name the file: {e}"
        );
    }

    #[test]
    fn multi_line_and_commented_tree_file() {
        let dir = std::env::temp_dir().join(format!("ip_tree2_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("t2.nwk");
        std::fs::write(&p, "# my tree\n((A,B),\n (C,D));\n").unwrap();
        let t = Tree::read(&p).unwrap();
        assert_eq!(t.n_tips(), 4);
    }

    #[test]
    fn postorder_has_children_before_parents() {
        let t = Tree::from_newick("(((A,B),C),(D,E));").unwrap();
        let po = t.postorder();
        let pos = |n: NodeId| po.iter().position(|&x| x == n).unwrap();
        for (p, c) in t.edges() {
            assert!(pos(c) < pos(p), "child {c} before parent {p}");
        }
        assert_eq!(*po.last().unwrap(), t.root());
        assert_eq!(po.len(), t.len());
    }
}
