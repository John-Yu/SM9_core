#![allow(dead_code)]
use core::ops::{Index, IndexMut};

/// A buffer to hold values of size 2 * N. This is mostly
/// a hack that's necessary until `generic_const_exprs` is stable.
#[derive(Copy, Clone)]
#[repr(C, align(8))]
pub struct MulBuffer<const N: usize> {
    pub b0: [u64; N],
    pub b1: [u64; N],
}
impl<const N: usize> MulBuffer<N> {
    const fn new(b0: [u64; N], b1: [u64; N]) -> Self {
        Self { b0, b1 }
    }

    pub const fn zeroed() -> Self {
        let b = [0u64; N];
        Self::new(b, b)
    }
    #[inline(always)]
    pub const fn get(&self, index: usize) -> &u64 {
        if index < N {
            &self.b0[index]
        } else {
            &self.b1[index - N]
        }
    }

    #[inline(always)]
    pub fn get_mut(&mut self, index: usize) -> &mut u64 {
        if index < N {
            &mut self.b0[index]
        } else {
            &mut self.b1[index - N]
        }
    }
}
impl<const N: usize> Index<usize> for MulBuffer<N> {
    type Output = u64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        self.get(index)
    }
}
impl<const N: usize> IndexMut<usize> for MulBuffer<N> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.get_mut(index)
    }
}
#[macro_export]
macro_rules! mac_with_carry {
    ($a:expr, $b:expr, $c:expr, &mut $carry:expr$(,)?) => {{
        let tmp = ($a as u128) + ($b as u128 * $c as u128) + ($carry as u128);
        $carry = (tmp >> 64) as u64;
        tmp as u64
    }};
}

#[macro_export]
macro_rules! adc {
    ($a:expr, $b:expr, &mut $carry:expr$(,)?) => {{
        let tmp = ($a as u128) + ($b as u128) + ($carry as u128);
        $carry = (tmp >> 64) as u64;
        tmp as u64
    }};
}
/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}
/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}
/// Calculate a + b * c, discarding the lower 64 bits of the result and setting
/// `carry` to the upper 64 bits.
#[inline(always)]
pub fn mac_discard(a: u64, b: u64, c: u64, carry: &mut u64) {
    let tmp = (a as u128) + (b as u128 * c as u128);
    *carry = (tmp >> 64) as u64;
}
/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}
