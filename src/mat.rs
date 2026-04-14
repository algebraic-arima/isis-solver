use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use bitvec::prelude::*;
use rand::Rng;
use std::ops::{Add, Div, Sub};

use crate::mod2mat::mod2matrix;

/// Matrix represented as a contiguous column-major buffer over Z_q.
/// data[c * row + r] means row r of column c.
#[derive(Debug, Clone)]
pub struct ColMatrixModQ {
    pub(crate) data: Vec<usize>,
    pub row: usize,
    pub(crate) col: usize,
    pub q: usize,
}

impl ColMatrixModQ {
    pub fn new(row: usize, q: usize) -> Self {
        if q == 0 {
            panic!("modulus q must be > 0");
        }
        Self {
            data: Vec::new(),
            row,
            col: 0,
            q,
        }
    }

    pub fn rand(row: usize, col: usize, q: usize) -> Self {
        if q == 0 {
            panic!("modulus q must be > 0");
        }
        let mut rng = rand::thread_rng();
        let mut data = vec![0_usize; row * col];
        for v in data.iter_mut() {
            *v = rng.gen_range(0..q);
        }
        Self { data, row, col, q }
    }

    pub fn append_col(&mut self, col_vec: VectorModQ) {
        if col_vec.data.len() != self.row {
            panic!("column length must equal row count");
        }
        self.data.reserve(self.row);
        self.data
            .extend(col_vec.data.into_iter().map(|x| x % self.q));
        self.col += 1;
    }

    pub fn trans_to_mod2mat(&self, ind: usize) -> mod2matrix {
        println!("trans_to_mod2mat: ind={}, row={}, col={}", ind, self.row, self.col);
        let mut out = mod2matrix::new(self.row);
        let n = self.row;
        for r in 0..n {
            for c in 0..(2 * n + 1) {
                out.set(r, c, self[(r, c + ind * (2 * n + 1))] % 2 == 1);
            }
        }
        out
    }

    pub fn linear_combine_01_bits(
        &self,
        selector: &BitVec<usize, Lsb0>,
        offset: usize,
    ) -> VectorModQ {
        if offset + selector.len() > self.col {
            panic!(
                "selector length exceeds: max {}, got {}",
                self.col,
                offset + selector.len()
            );
        }

        let mut out = vec![0_usize; self.row];
        for j in selector.iter_ones() {
            let start = (j + offset) * self.row;
            let col_slice = &self.data[start..start + self.row];
            for (out_i, &v) in out.iter_mut().zip(col_slice.iter()) {
                *out_i = (*out_i + v) % self.q;
            }
        }
        VectorModQ {
            data: out,
            row: self.row,
            q: self.q,
        }
    }
}

impl Index<(usize, usize)> for ColMatrixModQ {
    type Output = usize;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (r, c) = index;
        if r >= self.row || c >= self.col {
            panic!("index out of bounds");
        }
        &self.data[r + c * self.row]
    }
}

impl IndexMut<(usize, usize)> for ColMatrixModQ {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (r, c) = index;
        if r >= self.row || c >= self.col {
            panic!("index out of bounds");
        }
        &mut self.data[r + c * self.row]
    }
}

impl Display for ColMatrixModQ {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ColMatrixModQ ({}x{}, q={}):\n",
            self.row, self.col, self.q
        )?;
        for r in 0..self.row {
            write!(f, "[ ")?;
            for c in 0..self.col {
                write!(f, "{} ", self[(r, c)])?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}

impl Div<usize> for ColMatrixModQ {
    type Output = Self;

    fn div(self, rhs: usize) -> Self::Output {
        if rhs == 0 {
            panic!("division by zero");
        }
        assert!(
            self.q % rhs == 0,
            "division only supported when modulus is divisible by rhs"
        );
        let factor = self.q / rhs;
        let mut data = self.data;
        for v in data.iter_mut() {
            assert!(*v % rhs == 0, "matrix entry not divisible by rhs");
            *v /= rhs;
        }
        Self {
            data,
            row: self.row,
            col: self.col,
            q: factor,
        }
    }
}

#[derive(Debug, Clone)]
pub struct VectorModQ {
    data: Vec<usize>,
    row: usize,
    q: usize,
}

impl PartialEq for VectorModQ {
    fn eq(&self, other: &Self) -> bool {
        self.row == other.row && self.q == other.q && self.data == other.data
    }
}

impl VectorModQ {
    pub fn new(row: usize, q: usize) -> Self {
        if q == 0 {
            panic!("modulus q must be > 0");
        }
        Self {
            data: vec![0; row],
            row,
            q,
        }
    }

    pub fn from_bitvec(bitvec: &BitVec<usize, Lsb0>, q: usize) -> Self {
        if q == 0 {
            panic!("modulus q must be > 0");
        }
        let data = bitvec.iter().map(|b| if *b { 1 } else { 0 }).collect();
        Self {
            data,
            row: bitvec.len(),
            q,
        }
    }

    pub fn rand(row: usize, q: usize) -> Self {
        if q == 0 {
            panic!("modulus q must be > 0");
        }
        let mut rng = rand::thread_rng();
        Self {
            data: (0..row).map(|_| rng.gen_range(0..q)).collect(),
            row,
            q,
        }
    }

    pub fn rows(&self) -> usize {
        self.row
    }

    pub fn modulus(&self) -> usize {
        self.q
    }

    pub fn random_split_mod2(&self, m: usize) -> Vec<BitVec<usize, Lsb0>> {
        if m == 0 {
            panic!("m must be > 0");
        }

        let mut rng = rand::thread_rng();
        let n = self.data.len();

        let parity: BitVec<usize, Lsb0> = self.data.iter().map(|&v| v % 2 == 1).collect();

        let mut shares = vec![bitvec![usize, Lsb0; 0; n]; m];

        for j in 0..n {
            let mut acc = parity[j];
            for i in 0..(m - 1) {
                shares[i].set(j, rng.gen_bool(0.5));
                acc ^= shares[i][j];
            }
            shares[m - 1].set(j, acc);
        }
        shares
    }

    pub fn into_bitvec(&self) -> BitVec<usize, Lsb0> {
        assert!(self.q == 2, "vector modulus must be 2");
        self.data.iter().map(|&v| v % 2 == 1).collect()
    }
}

impl Display for VectorModQ {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[ ")?;
        for r in 0..self.row {
            write!(f, "{} ", self[r])?;
        }
        writeln!(f, "]")?;
        Ok(())
    }
}

impl Add for VectorModQ {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.row != rhs.row || self.q != rhs.q {
            panic!("vector dimensions or modulus do not match");
        }
        let data = self
            .data
            .into_iter()
            .zip(rhs.data.into_iter())
            .map(|(a, b)| (a + b) % self.q)
            .collect();
        Self {
            data,
            row: self.row,
            q: self.q,
        }
    }
}

impl Sub for VectorModQ {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.row != rhs.row || self.q != rhs.q {
            panic!("vector dimensions or modulus do not match");
        }
        let data = self
            .data
            .into_iter()
            .zip(rhs.data.into_iter())
            .map(|(a, b)| (a + self.q - b) % self.q)
            .collect();
        Self {
            data,
            row: self.row,
            q: self.q,
        }
    }
}

impl Div<usize> for VectorModQ {
    type Output = Self;

    fn div(self, rhs: usize) -> Self::Output {
        if rhs == 0 {
            panic!("division by zero");
        }
        assert!(
            self.q % rhs == 0,
            "division only supported when modulus is divisible by rhs"
        );
        let factor = self.q / rhs;
        let mut data = self.data;
        for i in 0..data.len() {
            assert!(data[i] % rhs == 0, "vector entry not divisible by rhs");
            data[i] /= rhs;
        }
        Self {
            data,
            row: self.row,
            q: factor,
        }
    }
}

impl Sub<&VectorModQ> for &VectorModQ {
    type Output = VectorModQ;

    fn sub(self, rhs: &VectorModQ) -> VectorModQ {
        if self.row != rhs.row || self.q != rhs.q {
            panic!("vector dimensions or modulus do not match");
        }
        let data = self
            .data
            .iter()
            .zip(rhs.data.iter())
            .map(|(&a, &b)| (a + self.q - b) % self.q)
            .collect();
        VectorModQ {
            data,
            row: self.row,
            q: self.q,
        }
    }
}

impl Div<usize> for &VectorModQ {
    type Output = VectorModQ;

    fn div(self, rhs: usize) -> VectorModQ {
        if rhs == 0 {
            panic!("division by zero");
        }
        assert!(
            self.q % rhs == 0,
            "division only supported when modulus is divisible by rhs"
        );
        let factor = self.q / rhs;
        let data: Vec<usize> = self
            .data
            .iter()
            .map(|&v| {
                assert!(v % rhs == 0, "vector entry not divisible by rhs");
                v / rhs
            })
            .collect();
        VectorModQ {
            data,
            row: self.row,
            q: factor,
        }
    }
}

impl Index<usize> for VectorModQ {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.row {
            panic!("index out of bounds");
        }
        &self.data[index]
    }
}

impl IndexMut<usize> for VectorModQ {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index >= self.row {
            panic!("index out of bounds");
        }
        &mut self.data[index]
    }
}
