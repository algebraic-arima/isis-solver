use bitvec::prelude::*;
use rand::Rng;
use std::fmt;

use crate::mat::ColMatrixModQ;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone)]
pub struct mod2matrix {
    pub data: BitVec<usize, Lsb0>, // row-major
    pub row: usize,
    pub col: usize,
}

impl mod2matrix {
    #[inline]
    fn index_of(&self, r: usize, c: usize) -> usize {
        r * self.col + c
    }

    fn swap_rows(&mut self, r1: usize, r2: usize) {
        if r1 == r2 {
            return;
        }
        for c in 0..self.col {
            let i1 = self.index_of(r1, c);
            let i2 = self.index_of(r2, c);
            self.data.swap(i1, i2);
        }
    }

    fn xor_row_from_col(&mut self, target: usize, source: usize, from_col: usize) {
        for c in from_col..self.col {
            let ti = self.index_of(target, c);
            let si = self.index_of(source, c);
            let next = self.data[ti] ^ self.data[si];
            self.data.set(ti, next);
        }
    }

    pub fn new(row: usize) -> Self {
        let col = row * 2 + 1;
        Self {
            data: bitvec![usize, Lsb0; 0; row * col],
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
        let mut data = bitvec![usize, Lsb0; 0; mat.row * mat.col];
        for r in 0..mat.row {
            for c in 0..mat.col {
                data.set(r * mat.col + c, mat.data[c * mat.row + r] % 2 == 1);
            }
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

        let mut data = bitvec![usize, Lsb0; 0; row * col];
        for (r, row_bits) in rows.into_iter().enumerate() {
            for c in 0..col {
                data.set(r * col + c, row_bits[c]);
            }
        }

        Self { data, row, col }
    }

    pub fn set(&mut self, r: usize, c: usize, value: bool) {
        if r >= self.row || c >= self.col {
            panic!(
                "Index out of bounds: ({}, {}) for matrix of size {}x{}",
                r, c, self.row, self.col
            );
        }
        let i = self.index_of(r, c);
        self.data.set(i, value);
    }

    pub fn get(&self, r: usize, c: usize) -> bool {
        if r >= self.row || c >= self.col {
            panic!("Index out of bounds");
        }
        self.data[self.index_of(r, c)]
    }

    pub fn append_row(&mut self, row: BitVec<usize, Lsb0>) {
        if row.len() != self.col {
            panic!("Row length does not match matrix width");
        }
        self.data.extend_from_bitslice(row.as_bitslice());
        self.row += 1;
    }

    pub fn rank(&self) -> usize {
        let mut mat = self.clone();
        let mut rank = 0;

        for col in 0..mat.col {
            if rank >= mat.row {
                break;
            }

            let pivot = (rank..mat.row).find(|&r| mat.get(r, col));
            if let Some(pivot_row) = pivot {
                mat.swap_rows(rank, pivot_row);

                for r in 0..mat.row {
                    if r != rank && mat.get(r, col) {
                        mat.xor_row_from_col(r, rank, col);
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

            if let Some(i) = (pivot_row..self.row).find(|&r| mat.get(r, col)) {
                mat.swap_rows(pivot_row, i);
                b.swap(pivot_row, i);

                let b_p = b[pivot_row];
                for r in (pivot_row + 1)..self.row {
                    if mat.get(r, col) {
                        mat.xor_row_from_col(r, pivot_row, col);
                        let br = b[r];
                        b.set(r, br ^ b_p);
                    }
                }

                pivot_pos[pivot_row] = Some(col);
                pivot_row += 1;
            }
        }
        assert!(pivot_row == self.row, "Matrix is not full row rank");

        let mut x = bitvec![usize, Lsb0; 0; self.col];
        for i in (0..pivot_row).rev() {
            if let Some(c) = pivot_pos[i] {
                // Σ (mat[i][j] * x[j]) for j > c
                let mut sum = false;
                for j in (c + 1)..self.col {
                    sum ^= mat.get(i, j) & x[j];
                }

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

            if let Some(i) = (pivot_row..self.row).find(|&r| mat.get(r, col)) {
                mat.swap_rows(pivot_row, i);
                b.swap(pivot_row, i);

                let b_p = b[pivot_row];
                for r in (pivot_row + 1)..self.row {
                    if mat.get(r, col) {
                        mat.xor_row_from_col(r, pivot_row, col);
                        let br = b[r];
                        b.set(r, br ^ b_p);
                    }
                }

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
                        sum ^= mat.get(i, j) & x[j];
                    } else {
                        assert!(j == c + 1);
                        x.set(j, false);
                        vis.set(j, true);
                        sum ^= mat.get(i, j) & x[j];
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
                        sum ^= mat.get(i, j) & y[j];
                    } else {
                        assert!(j == c + 1);
                        y.set(j, true);
                        vis.set(j, true);
                        sum ^= mat.get(i, j) & y[j];
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
                write!(f, "{} ", if self.get(r, c) { 1 } else { 0 })?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}
