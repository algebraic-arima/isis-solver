use bitvec::prelude::*;
use rand::Rng;
use std::fmt;

use crate::mat::ColMatrixModQ;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone)]
pub struct mod2matrix {
    pub data: Vec<BitVec<usize, Lsb0>>, // a row is a bitvec
    pub row: usize,
    pub col: usize,
}

impl mod2matrix {
    pub fn new(row: usize) -> Self {
        let col = row * 2 + 1;
        Self {
            data: (0..row).map(|_| bitvec![usize, Lsb0; 0; col]).collect(),
            row,
            col,
        }
    }

    pub fn from_mat_modq(mat: &ColMatrixModQ) -> Self {
        if mat.col != 2 * mat.row + 1 {
            panic!("Input matrix must be n x (2n+1)");
        }
        if mat.q != 2 {
            panic!("Input matrix modulus must be 2");
        }
        let mut data = Vec::with_capacity(mat.row);
        for r in 0..mat.row {
            let mut row_vec = bitvec![usize, Lsb0; 0; mat.col];
            for c in 0..mat.col {
                row_vec.set(c, mat.data[c * mat.row + r] % 2 == 1);
            }
            data.push(row_vec);
        }
        Self {
            data,
            row: mat.row,
            col: mat.col,
        }
    }

    pub fn from_rows(rows: Vec<BitVec<usize, Lsb0>>) -> Self {
        let row = rows.len();
        let col = rows.first().map_or(0, |r| r.len());

        if rows.iter().any(|r| r.len() != col) {
            panic!("All rows must have the same length");
        }

        Self {
            data: rows,
            row,
            col,
        }
    }

    pub fn set(&mut self, r: usize, c: usize, value: bool) {
        if r >= self.row || c >= self.col {
            panic!(
                "Index out of bounds: ({}, {}) for matrix of size {}x{}",
                r, c, self.row, self.col
            );
        }
        self.data[r].set(c, value);
    }

    pub fn get(&self, r: usize, c: usize) -> bool {
        if r >= self.row || c >= self.col {
            panic!("Index out of bounds");
        }
        self.data[r][c]
    }

    pub fn append_row(&mut self, row: BitVec<usize, Lsb0>) {
        if row.len() != self.col {
            panic!("Row length does not match matrix width");
        }
        self.data.push(row);
        self.row += 1;
    }

    pub fn rank(&self) -> usize {
        let mut mat = self.clone();
        let mut rank = 0;

        for col in 0..mat.col {
            if rank >= mat.row {
                break;
            }

            let pivot = (rank..mat.row).find(|&r| mat.data[r][col]);
            if let Some(pivot_row) = pivot {
                if pivot_row != rank {
                    mat.data.swap(rank, pivot_row);
                }

                let pivot_snapshot = mat.data[rank].clone();
                for r in 0..mat.row {
                    if r != rank && mat.data[r][col] {
                        for c in col..mat.col {
                            let next = mat.data[r][c] ^ pivot_snapshot[c];
                            mat.data[r].set(c, next);
                        }
                    }
                }

                rank += 1;
            }
        }

        rank
    }

    pub fn is_full_rank(&self) -> bool {
        self.rank() == self.row
    }

    fn basis_row(col: usize, index: usize) -> BitVec<usize, Lsb0> {
        let mut row = bitvec![usize, Lsb0; 0; col];
        if index < col {
            row.set(index, true);
        }
        row
    }

    pub fn solve_full_rank(&self, target: &BitVec<usize, Lsb0>) -> BitVec<usize, Lsb0> {
        if target.len() != self.row {
            panic!("Target vector length must match matrix row count");
        }

        let mut mat = self.clone();
        let mut b = target.clone();
        let mut pivot_row = 0;
        let mut pivot_pos = vec![None; self.row];

        for col in 0..self.col {
            if pivot_row >= self.row {
                break;
            }

            if let Some(i) = (pivot_row..self.row).find(|&r| mat.data[r][col]) {
                mat.data.swap(pivot_row, i);
                b.swap(pivot_row, i);

                let (upper, lower) = mat.data.split_at_mut(pivot_row + 1);
                let (b_upper, b_lower) = b.split_at_mut(pivot_row + 1);

                let r_p = &upper[pivot_row];
                let b_p = b_upper[pivot_row];

                lower
                    .iter_mut()
                    .zip(b_lower.iter_mut())
                    .for_each(|(row, mut row_b)| {
                        if row[col] {
                            *row ^= r_p;
                            *row_b ^= b_p;
                        }
                    });

                pivot_pos[pivot_row] = Some(col);
                pivot_row += 1;
            }
        }
        assert!(pivot_row == self.row, "Matrix is not full row rank");

        let mut x = bitvec![usize, Lsb0; 0; self.col];
        for i in (0..pivot_row).rev() {
            if let Some(c) = pivot_pos[i] {
                // Σ (mat[i][j] * x[j]) for j > c
                let sum = (mat.data[i][c + 1..].iter().by_vals())
                    .zip(x[c + 1..].iter().by_vals())
                    .fold(false, |acc, (m, xi)| acc ^ (m & xi));

                x.set(c, b[i] ^ sum);
            } else {
                panic!("Unexpected: pivot row without pivot column");
            }
        }

        x
    }

    pub fn procedure2(&self, y: &BitVec<usize, Lsb0>) -> BitVec<usize, Lsb0> {
        let n = self.row;
        if self.col != 2 * n + 1 {
            panic!("Matrix A must be n x (2n+1)");
        }

        let full_mat = self.extend_to_rank(2 * n + 1);

        let mut rng = rand::thread_rng();
        let mut u = bitvec![usize, Lsb0; 0; n + 1];
        for i in 0..(n + 1) {
            u.set(i, rng.gen_bool(0.5));
        }

        let mut combined_target = y.clone();
        combined_target.extend_from_bitslice(&u);

        full_mat.solve_full_rank(&combined_target)
    }

    pub fn extend_to_rank(&self, target_rank: usize) -> Self {
        let mut result = self.clone();
        let mut current_rank = result.rank();

        for i in 0..result.col {
            if current_rank >= target_rank {
                break;
            }

            let mut candidate = result.clone();
            candidate.append_row(Self::basis_row(result.col, i));
            let new_rank = candidate.rank();

            if new_rank > current_rank {
                result = candidate;
                current_rank = new_rank;
            }
        }
        result
    }

    pub fn solve_minus_one_rank(
        &self,
        target: &BitVec<usize, Lsb0>,
    ) -> (BitVec<usize, Lsb0>, BitVec<usize, Lsb0>) {
        if target.len() != self.row {
            panic!("Target vector length must match matrix row count");
        }

        let mut mat = self.clone();
        let mut b = target.clone();
        let mut pivot_row = 0;
        let mut pivot_pos = vec![None; self.row];

        for col in 0..self.col {
            if pivot_row >= self.row {
                break;
            }

            if let Some(i) = (pivot_row..self.row).find(|&r| mat.data[r][col]) {
                mat.data.swap(pivot_row, i);
                b.swap(pivot_row, i);

                let (upper, lower) = mat.data.split_at_mut(pivot_row + 1);
                let (b_upper, b_lower) = b.split_at_mut(pivot_row + 1);

                let r_p = &upper[pivot_row];
                let b_p = b_upper[pivot_row];

                lower
                    .iter_mut()
                    .zip(b_lower.iter_mut())
                    .for_each(|(row, mut row_b)| {
                        if row[col] {
                            *row ^= r_p;
                            *row_b ^= b_p;
                        }
                    });

                pivot_pos[pivot_row] = Some(col);
                pivot_row += 1;
            }
        }
        assert!(pivot_row == self.row, "Matrix is not full row rank");

        let mut x = bitvec![usize, Lsb0; 0; self.col];
        let mut vis = bitvec![usize, Lsb0; 0; self.col];
        for i in (0..pivot_row).rev() {
            if let Some(c) = pivot_pos[i] {
                // Σ (mat[i][j] * x[j]) for j > c
                let mut sum = false;
                for j in c + 1..self.col {
                    if vis[j] {
                        sum ^= mat.data[i][j] & x[j];
                    } else {
                        assert!(j == c + 1);
                        x.set(j, false);
                        vis.set(j, true);
                        sum ^= mat.data[i][j] & x[j];
                    }
                }
                x.set(c, b[i] ^ sum);
                vis.set(c, true);
            } else {
                panic!("Unexpected: pivot row without pivot column");
            }
        }

        let mut y = bitvec![usize, Lsb0; 0; self.col];
        let mut vis = bitvec![usize, Lsb0; 0; self.col];
        for i in (0..pivot_row).rev() {
            if let Some(c) = pivot_pos[i] {
                let mut sum = false;
                for j in c + 1..self.col {
                    if vis[j] {
                        sum ^= mat.data[i][j] & y[j];
                    } else {
                        assert!(j == c + 1);
                        y.set(j, true);
                        vis.set(j, true);
                        sum ^= mat.data[i][j] & y[j];
                    }
                }
                y.set(c, b[i] ^ sum);
                vis.set(c, true);
            } else {
                panic!("Unexpected: pivot row without pivot column");
            }
        }

        (x, y)
    }

    pub fn procedure1(
        &self,
        y: &BitVec<usize, Lsb0>,
    ) -> (BitVec<usize, Lsb0>, BitVec<usize, Lsb0>) {
        let n = self.row;
        if self.col != 2 * n + 1 {
            panic!("Matrix A must be n x (2n+1)");
        }

        let full_mat = self.extend_to_rank(2 * n);

        let mut rng = rand::thread_rng();
        let mut u = bitvec![usize, Lsb0; 0; n];
        for i in 0..n {
            u.set(i, rng.gen_bool(0.5));
        }

        let mut combined_target = y.clone();
        combined_target.extend_from_bitslice(&u);

        full_mat.solve_minus_one_rank(&combined_target)
    }
}

impl fmt::Display for mod2matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Matrix over F_2 ({}x{}):", self.row, self.col)?;
        for r in 0..self.row {
            write!(f, "[ ")?;
            for c in 0..self.col {
                write!(f, "{} ", if self.data[r][c] { 1 } else { 0 })?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}
