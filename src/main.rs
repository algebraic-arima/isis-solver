use bitvec::prelude::*;
use isis_solver::{mat::*, mod2mat::*};

const N: usize = 8;
const L: usize = 3;
const B: usize = N * 2 + 1;
const Q: usize = 1 << L;
const M: usize = B.pow(L as u32);

fn isis_solve(
    l: usize,
    m: usize,
    a: &ColMatrixModQ,
    y: &VectorModQ,
) -> BitVec<usize, Lsb0> {
    if l == 1 && m == B {
        let a = mod2matrix::from_mat_modq(a);
        let y = y.into_bitvec();
        return a.procedure2(&y);
    }
    let mm = m / B;
    let yvec = y.random_split_mod2(mm);
    let mut x_all: BitVec<usize, Lsb0> = BitVec::new();
    let mut a_new = ColMatrixModQ::new(a.row, a.q);
    let mut z_list = Vec::with_capacity(mm);
    
    for i in 0..mm {
        let target = &yvec[i];
        let mod2block = a.trans_to_mod2mat(i);
        if !mod2block.is_full_rank() {
            panic!("Unexpected: mod2 block is not full rank");
        }
        let (x1, x2) = mod2block.procedure1(target);

        z_list.push((x1.clone(), x2.clone()));
        let ax1 = a.linear_combine_01_bits(&x1, i * B);
        let ax2 = a.linear_combine_01_bits(&x2, i * B);
        let az = ax2 - ax1;
        a_new.append_col(az);
        x_all.extend(x1);
    }
    
    let dy = a.linear_combine_01_bits(&x_all, 0);
    let y_new = (y - &dy) / 2;
    let a_new = a_new / 2;
    let xx = isis_solve(l - 1, mm, &a_new, &y_new);
    
    let mut ans = BitVec::with_capacity(xx.len() * B);
    for i in 0..xx.len() {
        if xx[i] {
            ans.extend(z_list[i].1.iter());
        } else {
            ans.extend(z_list[i].0.iter());
        }
    }
    ans
}

fn main() {
    let a = ColMatrixModQ::rand(N, M, Q);
    let y = VectorModQ::rand(N, Q);
    let xf = isis_solve(L, M, &a, &y);
    let check = a.linear_combine_01_bits(&xf, 0);
    // println!("Matrix A:\n{}", a);
    println!("Vector y:\n{}", y);
    // println!("ans: {}", xf);
    println!("check: {}", check);
    assert!(check == y);
}
