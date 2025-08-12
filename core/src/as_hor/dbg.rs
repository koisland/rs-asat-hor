/// Based on Ben Langmead's implementation
/// * https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_deBruijn.ipynb#scrollTo=E3sblU0B8W0n
use std::{
    collections::{BinaryHeap, HashMap, HashSet},
    fmt::{Debug, Display},
    hash::Hash,
};

use eyre::bail;

/// Edge mapping left node id to right node ids
type Edges = HashMap<usize, Vec<usize>>;
/// Mapping of node id to type.
type NodeIDs<'a, T> = HashMap<usize, &'a [T]>;
/// Mapping of type to node metadata.
type Nodes<'a, T> = HashMap<&'a [T], Node<'a, T>>;

#[derive(Debug, Clone, PartialEq, Eq)]
struct Node<'a, T> {
    id: usize,
    elem: &'a [T],
    nin: usize,
    nout: usize,
}

impl<'a, T: Hash> Hash for Node<'a, T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.elem.hash(state);
    }
}

impl<'a, T: Debug> Display for Node<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.elem)
    }
}

impl<'a, T> Node<'a, T> {
    fn new(elem: &'a [T], id: usize) -> Self {
        Node {
            id,
            elem,
            nin: 0,
            nout: 0,
        }
    }

    fn is_semi_balanced(&self) -> bool {
        self.nin.abs_diff(self.nout) == 1
    }

    fn is_balanced(&self) -> bool {
        self.nin == self.nout
    }
}

#[derive(Debug)]
pub struct Dbg<'a, T: PartialEq + Eq + Hash, const N: usize> {
    nodes: Nodes<'a, T>,
    node_ids: NodeIDs<'a, T>,
    node_counts: HashMap<&'a [T], usize>,
    edges: Edges,
    nsemi: usize,
    nbal: usize,
    nneither: usize,
    head: Option<usize>,
    tail: Option<usize>,
}
#[derive(Debug, PartialEq, Eq)]
struct NodeCount<'a, T>(&'a [T], usize);

impl<'a, T> Ord for NodeCount<'a, T>
where
    T: PartialEq + Eq,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.1.cmp(&other.1)
    }
}
impl<'a, T> PartialOrd for NodeCount<'a, T>
where
    T: PartialEq + Eq,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.1.cmp(&other.1))
    }
}

impl<'a, T: PartialEq + Eq + Hash + Debug, const N: usize> Dbg<'a, T, N> {
    fn chop(elems: &[T]) -> impl Iterator<Item = (&[T], &[T], &[T])> {
        (0..elems.len() - (N - 1))
            .map(move |i| (&elems[i..i + N], &elems[i..i + N - 1], &elems[i + 1..i + N]))
    }

    pub fn new(elems: &'a [T]) -> Self {
        let mut nodes: Nodes<T> = HashMap::new();
        let mut node_ids: NodeIDs<T> = HashMap::new();
        let mut edges: Edges = HashMap::new();
        let mut n_node_ids = 0;
        let mut node_counts = HashMap::new();

        for (elem, elem_1l, elem_1r) in Dbg::<T, N>::chop(elems) {
            node_counts
                .entry(elem_1l)
                .and_modify(|x| *x += 1)
                .or_insert(1);

            let node_id_l = {
                let node_l = nodes.entry(elem_1l).or_insert_with(|| {
                    n_node_ids += 1;
                    Node::new(elem_1l, n_node_ids)
                });
                node_ids.insert(node_l.id, elem_1l);
                node_l.nout += 1;
                node_l.id
            };
            let node_id_r = {
                let node_r = nodes.entry(elem_1r).or_insert_with(|| {
                    n_node_ids += 1;
                    Node::new(elem_1r, n_node_ids)
                });
                node_ids.insert(node_r.id, elem_1r);
                node_r.nin += 1;
                node_r.id
            };

            edges
                .entry(node_id_l)
                .or_default()
                .push(node_id_r);
        }

        let mut nsemi = 0;
        let mut nbal = 0;
        let mut nneither = 0;
        let mut head = None;
        let mut tail = None;

        for node in nodes.values() {
            if node.is_balanced() {
                nbal += 1
            } else if node.is_semi_balanced() {
                if node.nin == node.nout + 1 {
                    tail = Some(node.id)
                }
                if Some(node.nin) == node.nout.checked_sub(1) {
                    head = Some(node.id)
                }
                nsemi += 1
            } else {
                nneither += 1
            }
        }
        Self {
            nodes,
            edges,
            node_ids,
            node_counts,
            nsemi,
            nbal,
            nneither,
            head,
            tail,
        }
    }

    fn nnodes(&self) -> usize {
        self.nodes.len()
    }

    fn nedges(&self) -> usize {
        self.edges.len()
    }

    /// Return true iff graph has Eulerian walk (visits every node once).
    fn has_eulerian_walk(&self) -> bool {
        self.nneither == 0 && self.nsemi == 2
    }

    /// Return true iff graph has Eulerian cycle (visits every node once and returns to same node).
    fn has_eulerian_cycle(&self) -> bool {
        self.nneither == 0 && self.nsemi == 0
    }

    /// Return true iff graph has Eulerian walk or cycle
    fn is_eulerian(&self) -> bool {
        // technically, if it has an Eulerian walk
        self.has_eulerian_walk() || self.has_eulerian_cycle()
    }

    fn find_cycles(&self) -> eyre::Result<Vec<Vec<usize>>> {
        // Iterate from largest node and greedily take next.
        // Similar to SRF's algo but we're operating at the monomer scale.
        let mut node_counts =
            BinaryHeap::from_iter(self.node_counts.iter().map(|nc| NodeCount(nc.0, *nc.1)));
        const MIN_OCC: usize = 2;

        let mut searches = vec![];
        while let Some(mut curr_node) = node_counts
            .pop()
            .and_then(|node| {
                node.1
                    .gt(&MIN_OCC)
                    .then(|| self.nodes.get(node.0))
                    .flatten()
            })
            .map(|node| node.id)
        {
            let starting_node = curr_node;
            let mut search = vec![];
            let mut is_cycle = false;
            let mut traveled_nodes = HashSet::new();

            while let Some(next_node) = self.edges
                .get(&curr_node)
                .and_then(|choices|
                    // Choose the largest occuring node at a bifurcation in graph.
                    choices.iter()
                    .filter(|node| **node == starting_node || !traveled_nodes.contains(node))
                    .max_by(|node_a, node_b|
                        self.node_counts[self.node_ids[node_a]].cmp(&self.node_counts[self.node_ids[node_b]])
                    )
                ) {
                    // Track traveled nodes.
                    traveled_nodes.insert(next_node);

                    // Hit end of cycle.
                    if *next_node == starting_node {
                        search.push(*next_node);
                        is_cycle = true;
                        break;
                    }
                    curr_node = *next_node;
                    search.push(*next_node);
                }

            if is_cycle {
                searches.push(search);
            }
        }

        Ok(searches)
    }

    // fn eularian_walk_or_cycle(&self) -> Option<Vec<[&T; N]>> {
    fn eularian_walk_or_cycle(&self) -> Option<Vec<&[T]>> {
        assert!(self.is_eulerian(), "Not eularian HOR.");

        let mut graph = self.edges.clone();
        // Ensure end to walk by adding head to tail?
        if self.has_eulerian_walk() {
            let (Some(tail), Some(head)) = (self.tail, self.head) else {
                return None;
            };
            graph.entry(tail).or_default().push(head);
        }

        // Graph now has eularian cycle.
        let mut search = vec![];

        // Get starting node.
        let mut curr_node = graph.keys().next().cloned()?;

        while let Some(next_node) = graph.get_mut(&curr_node).and_then(|choices| choices.pop()) {
            curr_node = next_node;
            search.push(next_node);
        }

        // Reverse and take all but last node.
        search.reverse();
        search.pop();

        // Adjust node list so that it starts at head and ends at tail
        let search = if self.has_eulerian_walk() {
            let head_idx = search.iter().position(|idx| Some(*idx) == self.head)?;
            search
                .get(head_idx..)
                .unwrap()
                .iter()
                .chain(search.get(..head_idx).unwrap())
                .cloned()
                .collect()
        } else {
            search
        };

        Some(
            search
                .into_iter()
                .filter_map(|nid| self.node_ids.get(&nid))
                .cloned()
                .collect(),
        )
    }
}

#[cfg(test)]
mod test {
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufWriter, Write},
    };

    use dot_structures::*;
    use graphviz_rust::printer::{DotPrinter, PrinterContext};
    use itertools::Itertools;

    use crate::{as_hor::dbg::Dbg, HOR};

    fn hor_repeating() -> Vec<String> {
        const HOR: &str = "S2C4H1L.5-14_8-9_3-14_8-9_3-14_8-14_8-9_3-14_8-9_3-14_8-9_3-14_8-9_3-14_8-10_4-14_8-9_3-14_8-14_8-9_3-14_8-9_3-14_8-9_3-14_8-9_3-14_8-9_3-19";
        HOR::new(HOR)
            .unwrap()
            .monomers()
            .iter()
            .map(|m| {
                let mut mon = m.monomers.iter().join("/");
                mon.push('_');
                mon
            })
            .collect_vec()
    }

    // #[test]
    // fn test_eularian_walk_cycle() {
    //     let monomers = hors_cen();
    //     // let monomers = hor_repeating();
    //     let dbg = Dbg::<String, 7>::new(&monomers);
    //     let res = dbg.eularian_walk_or_cycle().unwrap();
    //     println!(
    //         "{:?}",
    //         res.iter().map(|e| e.join("")).collect_vec()
    //     )
    // }

    #[test]
    fn test_print_dot() {
        let monomers = hor_repeating();
        let dbg = Dbg::<String, 10>::new(&monomers);

        // let mut writer = BufWriter::new(File::create("out.dot").unwrap());
        let mut dot_dbg = Graph::DiGraph {
            id: Id::Plain("DBG".to_owned()),
            strict: true,
            stmts: vec![],
        };
        for node_id in dbg.edges.keys() {
            let node = dbg.node_ids[node_id];
            if dbg.node_counts.get(node) < Some(&20) {
                continue;
            }
            let node_str = format!("\"{}\"", node.join("").trim_end_matches('_'));
            dot_dbg.add_stmt(Stmt::Node(Node::new(
                NodeId(Id::Plain(node_str.clone()), None),
                vec![],
            )));
        }
        for (node_id, dsts) in dbg.edges.iter() {
            let node = dbg.node_ids[node_id];
            let node_str = format!("\"{}\"", node.join("").trim_end_matches('_'));
            let mut wt_map: HashMap<usize, usize> = HashMap::new();
            for dst in dsts {
                wt_map.entry(*dst).and_modify(|c| *c += 1).or_insert(1);
            }
            for (dst_node_id, dst_node_wt) in wt_map.iter().filter(|(a, b)| **b > 20) {
                let dst_node = dbg.node_ids[dst_node_id];
                let dst_node_str = format!("\"{}\"", dst_node.join("").trim_end_matches('_'));
                dot_dbg.add_stmt(Stmt::Edge(Edge {
                    ty: EdgeTy::Pair(
                        Vertex::N(NodeId(Id::Plain(node_str.clone()), None)),
                        Vertex::N(NodeId(Id::Plain(dst_node_str), None)),
                    ),
                    attributes: vec![Attribute(
                        Id::Plain(String::from("label")),
                        Id::Plain(format!("{dst_node_wt}")),
                    )],
                }));
            }
        }

        let mut ctx = PrinterContext::default();
        ctx.with_indent_step(4);
        let dot_dbg_str = dot_dbg.print(&mut ctx);
        // writeln!(&mut writer, "{}", dot_dbg_str).unwrap();
    }
}
