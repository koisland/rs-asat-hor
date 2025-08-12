use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Write},
};

use itertools::Itertools;
use rust_lapper::{Interval, Lapper};
use suffix::SuffixTable;

fn hors_cen() -> (String, Vec<(usize, usize)>, HashMap<char, String>) {
    let reader = BufReader::new(File::open("test/final_decomposition.bed").unwrap());

    let mut monomer_seq = String::new();
    let mut monomer_coords = vec![];
    let mut monomer_counts: HashMap<String, usize> = HashMap::new();
    let mut char_monomer_map = HashMap::new();
    let mut monomer_char_map: HashMap<String, char> = HashMap::new();

    for line in reader.lines().map_while(Result::ok) {
        let Some((chrom, st, end, mon, score, strand, tst, tend, _)) =
            line.split('\t').collect_tuple()
        else {
            continue;
        };
        let score = score.parse::<f32>().unwrap();
        if score < 70.0 {
            continue;
        }

        let monomer_char = if let Some(mon_cnt) = monomer_counts.get_mut(mon) {
            *mon_cnt += 1;
            monomer_char_map[mon]
        } else {
            monomer_counts.insert(mon.to_owned(), 1);
            // Store monomer as char starting from char(1)
            let monomer_char = char::from_u32(monomer_counts.len() as u32)
                .expect("Overflowed. Too many monomers.");
            char_monomer_map.insert(monomer_char, mon.to_owned());
            monomer_char_map.insert(mon.to_owned(), monomer_char);
            monomer_char
        };
        monomer_seq.push(monomer_char);
        monomer_coords.push((st.parse::<usize>().unwrap(), end.parse::<usize>().unwrap()));
    }
    (monomer_seq, monomer_coords, char_monomer_map)
}

fn rle_counts(seq: &str) -> Vec<(char, i32)> {
    seq.chars().fold(vec![], |mut acc: Vec<(char, i32)>, x| {
        if let Some(mon) = acc.last_mut().filter(|mon| mon.0 == x) {
            mon.1 += 1
        } else {
            acc.push((x, 1));
        }
        acc
    })
}

fn format_rle_counts(counts: &[(char, i32)], char_to_str: &HashMap<char, String>) -> String {
    counts
        .iter()
        .map(|(c, cnt)| {
            if *cnt > 1 {
                format!("{}.{cnt}-", char_to_str[c])
            } else {
                format!("{}-", char_to_str[c])
            }
        })
        .join("")
}

#[test]
pub fn detect_repeats() {
    let (seq, seq_coords, char_to_mon) = hors_cen();
    // let mut writer = std::io::BufWriter::new(File::create("out.txt").unwrap());

    // Construct the suffix table and longest common prefix array.
    let sfx_tbl = SuffixTable::new(&seq);
    let lcp_arr = sfx_tbl.lcp_lens();

    let mut all_repeats: Vec<Interval<usize, String>> = vec![];
    for (idx_sfx, sfx_length) in lcp_arr.into_iter().enumerate() {
        let repeat = &sfx_tbl.suffix(idx_sfx)[0..sfx_length as usize];

        let positions = sfx_tbl.positions(repeat);
        let mut positions_iter = positions.iter().sorted().peekable();
        let mut total_length = 0;
        let mut differences = vec![];
        let mut valid_positions = vec![];
        loop {
            let Some(pos) = positions_iter.next() else {
                break;
            };
            let Some(next_pos) = positions_iter.peek() else {
                total_length += sfx_length;
                break;
            };
            let diff = *next_pos - pos;
            match diff.cmp(&sfx_length) {
                // Some overlap if diff between two adjacent positions less than the largest sfx length.
                Ordering::Less => {
                    differences.push(diff);
                    valid_positions.push(pos);
                    total_length += diff;
                }
                Ordering::Equal => {
                    differences.push(diff);
                    valid_positions.push(pos);
                    total_length += sfx_length
                }
                // But if diff is larger, indicates suffixes are not adjacent and should be ignored in total length calculation.
                Ordering::Greater => {
                    continue;
                }
            }
        }
        // Not a repeat. Single unit.
        if total_length == sfx_length {
            continue;
        }

        let Some(repeat_diff_mode) = differences
            .into_iter()
            .counts()
            .into_iter()
            .max_by(|a, b| a.1.cmp(&b.1))
            .map(|m| m.0)
        else {
            continue;
        };

        if let Some(smallest_repeat) = repeat.get(0..repeat_diff_mode as usize) {
            let smallest_repeat_counts = rle_counts(smallest_repeat);
            let mut min_rpt = format_rle_counts(&smallest_repeat_counts, &char_to_mon);
            min_rpt.pop();

            // let repeat_counts = rle_counts(repeat);
            // let mut conv_rpt = format_rle_counts(&repeat_counts, &char_to_mon);
            // conv_rpt.pop();

            all_repeats.extend(valid_positions.iter().map(|p| {
                let itv = seq_coords[**p as usize];
                Interval {
                    start: itv.0,
                    stop: itv.1 + 1,
                    val: min_rpt.to_owned(),
                }
            }));
        }
    }
    let mut all_repeats_itree = Lapper::new(all_repeats);
    all_repeats_itree.merge_overlaps();

    let chrom = "mPanPan1_chr17_chr19_pat_hsa17:38408651-45141749";
    for itv in all_repeats_itree.iter() {
        // writeln!(&mut writer, "{chrom}\t{}\t{}\t{}\t{}", itv.start, itv.stop, itv.val, itv.stop-itv.start).unwrap();
    }
    // Some(())
}
