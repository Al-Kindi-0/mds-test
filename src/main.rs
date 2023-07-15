use std::{collections::HashMap, fs};

use itertools::Itertools;
use ndarray::{array, Array2};

const MODULUS: u32 = 2_u32.pow(31) - 1;

// Reduced set of minors generation
fn compositions(k: usize, n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![];
    }
    if k == 1 {
        return vec![vec![n]];
    }

    let mut comp = vec![];
    for i in 1..=n {
        for t in compositions(k - 1, n - i) {
            let mut new_comp = vec![i];
            new_comp.extend(t);
            comp.push(new_comp);
        }
    }
    comp
}

fn left_shift(tup: &Vec<usize>, n: usize) -> Vec<usize> {
    if tup.is_empty() || n == 0 {
        return tup.clone();
    }
    let n = n % tup.len();
    let mut left = tup[n..].to_vec();
    left.extend(tup[..n].to_vec());
    left
}

fn generate_d(u: usize, n: usize) -> Vec<Vec<usize>> {
    let mut d = compositions(u, n);
    let mut d_shifts = Vec::new();

    let mut out = vec![];
    while !d.is_empty() {
        let t = d.remove(0);
        d_shifts.clear();
        for k in 0..u {
            d_shifts.push(left_shift(&t.clone(), k));
        }
        for dup in &d_shifts {
            if d.contains(dup) {
                d.retain(|x| x != dup);
            }
        }
        out.push(t);
    }

    out
}

fn accumu(lis: &[usize]) -> Vec<usize> {
    let mut total = 0;
    lis.iter()
        .map(|x| {
            total += x;
            total
        })
        .collect()
}

fn generate_i(d: &Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut result = vec![];
    for u in d {
        let mut tmp: Vec<usize> = accumu(&[0].to_vec().iter().chain(u.iter()).cloned().collect::<Vec<usize>>());
        tmp.pop();
        result.push(tmp);
    }
    result
}

fn find_subsets(s: Vec<usize>, n: usize) -> Vec<Vec<usize>> {
    s.into_iter().combinations(n).collect()
}

fn generate_j(u: usize, n: usize) -> Vec<Vec<usize>> {
    let j: Vec<usize> = (0..n as usize).collect();
    find_subsets(j, u)
}

fn circular_shifts_mod(i: &Vec<usize>, n: usize) -> Vec<Vec<usize>> {
    let mut result = vec![];
    let _u = i.len();
    for k in 0..n {
        let mut tmp = vec![];
        for ii in i {
            tmp.push((k + ii) % n);
        }
        tmp.sort();
        result.push(tmp);
    }
    result
}

fn minimal_set_of_minors_incremental(dim: usize) {
    for u in 2..dim {
        println!("Generating minors of size {:?} for dimension {:?}", u, dim);
        let d = generate_d(u, dim);
        let i_u = generate_i(&d);
        let mut j_u = generate_j(u, dim);
        let mut j_p = vec![];
        let mut m_sub = vec![];
        for i in i_u.iter() {
            j_u = j_u.iter().filter(|&x| !j_p.contains(x)).cloned().collect();
            for j in &j_u {
                m_sub.push((i.clone(), j.clone()));
            }
            for tmp in circular_shifts_mod(&i, dim) {
                j_p.push(tmp);
            }
        }
        let strng = serde_json::to_string(&m_sub).unwrap();
        fs::write(format!("./minors_set/{dim}/{u}.txt"), strng).expect("Unable to write file");
    }
    let mut all_indices = Vec::new();
    let mut all_indices_ = Vec::new();
    for i in 0..dim {
        all_indices.push(i);
        all_indices_.push(i);
    }

    let all_min = (
        Vec::from(all_indices.clone()),
        Vec::from(all_indices.clone()),
    );
    let strng = serde_json::to_string(&vec![all_min]).unwrap();
    fs::write(format!("./minors_set/{dim}/{dim}.txt"), strng).expect("Unable to write file");
}

// MDS test
fn read_list_minors(dim: usize, k: usize) -> Vec<(Vec<usize>, Vec<usize>)> {
    let strng =
        fs::read_to_string(format!("./minors_set/{dim}/{k}.txt")).expect("Unable to read file");
    serde_json::from_str(&strng).unwrap()
}
pub fn test_mds_incremental_step(
    a: &Array2<u32>,
    dim: usize,
    step: usize,
    current_minors: &HashMap<(Vec<usize>, Vec<usize>), u32>,
) -> HashMap<(Vec<usize>, Vec<usize>), u32> {
    // read the minimal effective set of minors at step `step` from disk
    let m: Vec<(Vec<usize>, Vec<usize>)> = read_list_minors(dim, step);

    let p: u32 = MODULUS;
    let mut new_all_minors = HashMap::new();
    for m_item in m.iter() {
        let (rows, cols) = m_item;
        let last_row = rows.len() - 1;
        let i = rows[last_row];

        let mut det = 0;

        for c in 0..cols.len() {
            let mut tmp_col = cols.clone();
            tmp_col.remove(c);

            let cofactor;
            if c % 2 == 1 {
                cofactor = ((p - a[(i, cols[c])]) as u64)
                    * (current_minors[&(rows[0..last_row].to_vec(), tmp_col)] as u64);
            } else {
                cofactor = (a[(i, cols[c])] as u64)
                    * (current_minors[&(rows[0..last_row].to_vec(), tmp_col)] as u64);
            }

            det = ((det as u64 + (cofactor % p as u64)) % p as u64) as u32;
        }

        if det == 0 {
            println!("The failing minor is {:?}", m_item);
            panic!("Non MDS: Minor is zero")
        }
        new_all_minors.insert((rows.clone(), cols.clone()), det);
    }
    new_all_minors
}

pub fn test_mds_incremental(a: &Array2<u32>, dim: usize) -> bool {
    let mut current_minors = HashMap::new();
    for i in 0..dim {
        for j in 0..dim {
            if a[(i, j)] == 0 {
                return false;
            }
            current_minors.insert((vec![i], vec![j]), a[(i, j)]);
        }
    }
    let _p: u32 = 2_u32.pow(31) - 1;

    for i in 2..=dim {
        println!("Testing {:?} minors", i);
        current_minors = test_mds_incremental_step(&a, dim, i, &current_minors);
    }

    true
}

fn main() {
    let dim = 12;
    let a = array![
        [7, 8, 21, 22, 6, 7, 9, 10, 13, 26, 8, 23],
        [23, 7, 8, 21, 22, 6, 7, 9, 10, 13, 26, 8],
        [8, 23, 7, 8, 21, 22, 6, 7, 9, 10, 13, 26],
        [26, 8, 23, 7, 8, 21, 22, 6, 7, 9, 10, 13],
        [13, 26, 8, 23, 7, 8, 21, 22, 6, 7, 9, 10],
        [10, 13, 26, 8, 23, 7, 8, 21, 22, 6, 7, 9],
        [9, 10, 13, 26, 8, 23, 7, 8, 21, 22, 6, 7],
        [7, 9, 10, 13, 26, 8, 23, 7, 8, 21, 22, 6],
        [6, 7, 9, 10, 13, 26, 8, 23, 7, 8, 21, 22],
        [22, 6, 7, 9, 10, 13, 26, 8, 23, 7, 8, 21],
        [21, 22, 6, 7, 9, 10, 13, 26, 8, 23, 7, 8],
        [8, 21, 22, 6, 7, 9, 10, 13, 26, 8, 23, 7]
    ];

    minimal_set_of_minors_incremental(dim);
    test_mds_incremental(&a, dim);
}
