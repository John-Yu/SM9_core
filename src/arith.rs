#![allow(dead_code)]
/// Divide by two
#[inline]
pub fn div2(a: &mut [u128; 2]) {
    let tmp = a[1] << 127;
    a[1] >>= 1;
    a[0] >>= 1;
    a[0] |= tmp;
}
/// Multiply by two, return true if has carrry
#[inline]
pub fn mul2(a: &mut [u128; 2]) -> bool {
    let tmp = a[0] >> 127;
    let carry = a[1] >> 127;
    a[0] <<= 1;
    a[1] <<= 1;
    a[1] |= tmp;
    carry != 0
}

#[inline(always)]
pub fn split_u128(i: u128) -> (u128, u128) {
    (i >> 64, i & 0xFFFFFFFFFFFFFFFF)
}

#[inline(always)]
pub fn combine_u128(hi: u128, lo: u128) -> u128 {
    (hi << 64) | lo
}

#[inline]
pub fn adc(a: u128, b: u128, carry: &mut u128) -> u128 {
    let (a1, a0) = split_u128(a);
    let (b1, b0) = split_u128(b);
    let (c, r0) = split_u128(a0 + b0 + *carry);
    let (c, r1) = split_u128(a1 + b1 + c);
    *carry = c;

    combine_u128(r1, r0)
}
#[inline]
pub fn add_nocarry(a: &mut [u128; 2], b: &[u128; 2]) {
    let mut carry;

    (a[0], carry) = a[0].overflowing_add(b[0]);
    (a[1], carry) = carrying_add(a[1], b[1], carry);

    debug_assert!(!carry);
}
// add , return true if has carry
#[inline]
pub fn add_carry(a: &mut [u128; 2], b: &[u128; 2]) -> bool {
    let mut carry;

    (a[0], carry) = a[0].overflowing_add(b[0]);
    (a[1], carry) = carrying_add(a[1], b[1], carry);

    carry
}
#[inline]
pub fn sub_noborrow(a: &mut [u128; 2], b: &[u128; 2]) {
    let mut borrow;

    (a[0], borrow) = a[0].overflowing_sub(b[0]);
    (a[1], borrow) = borrowing_sub(a[1], b[1], borrow);

    debug_assert!(!borrow);
}
// sub , return true if has borrow
#[inline]
pub fn sub_borrow(a: &mut [u128; 2], b: &[u128; 2]) -> bool {
    let mut borrow;

    (a[0], borrow) = a[0].overflowing_sub(b[0]);
    (a[1], borrow) = borrowing_sub(a[1], b[1], borrow);

    borrow
}
/// this = this * by * R^-1 (mod modulus), where R = 2^256
// The Montgomery reduction here is based on Algorithm 14.32 in
// Handbook of Applied Cryptography
// <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.
#[inline]
pub fn mul_reduce(this: &mut [u128; 2], by: &[u128; 2], modulus: &[u128; 2], inv: u128) -> bool {
    let mut res = mul(this, by);
    let carry = reduce(&mut res, modulus, inv);
    this.copy_from_slice(&res[2..]);

    carry
}
/// this = this^2 * R^-1 (mod modulus), where R = 2^256
#[inline]
pub fn square_reduce(this: &mut [u128; 2], modulus: &[u128; 2], inv: u128) -> bool {
    let mut res = square(this);
    let carry = reduce(&mut res, modulus, inv);
    this.copy_from_slice(&res[2..]);

    carry
}
/// Calculates the “full multiplication” lhs * rhs + carry without the possibility to overflow.
/// This returns the low-order (wrapping) bits and the high-order (overflow) bits of the result as two separate values, in that order.
/// u128::carrying_mul is not existed, so use this version.
// come from crate awint.
#[inline]
pub(crate) fn carrying_mul(lhs: u128, rhs: u128, carry: u128) -> (u128, u128) {
    let (sum0, sum1) = widening_mul(lhs, rhs);
    let (sum0, carry2) = sum0.overflowing_add(carry);
    let sum1 = sum1.wrapping_add(carry2 as u128);

    (sum0, sum1)
}
/// Calculates the complete product lhs * rhs without the possibility to overflow.
/// This returns the low-order (wrapping) bits and the high-order (overflow) bits of the result as two separate values, in that order.
/// If you also need to add a carry to the wide result, then you want carrying_mul instead.
/// u128::widening_mul is not existed, so use this version.
#[inline]
pub(crate) fn widening_mul(lhs: u128, rhs: u128) -> (u128, u128) {
    //                       [rhs_hi]  [rhs_lo]
    //                       [lhs_hi]  [lhs_lo]
    //                     X___________________
    //                       [------tmp0------]
    //             [------tmp1------]
    //             [------tmp2------]
    //     [------tmp3------]
    // +_______________________________________
    //                       [------sum0------]
    //     [------sum1------]

    let (lhs_hi, lhs_lo) = split_u128(lhs);
    let (rhs_hi, rhs_lo) = split_u128(rhs);
    let tmp0 = lhs_lo.wrapping_mul(rhs_lo);
    let tmp1 = lhs_lo.wrapping_mul(rhs_hi);
    let tmp2 = lhs_hi.wrapping_mul(rhs_lo);
    let tmp3 = lhs_hi.wrapping_mul(rhs_hi);
    // tmp1 and tmp2 straddle the boundary. We have to handle three carries
    let (sum0, carry0) = tmp0.overflowing_add(tmp1.wrapping_shl(64));
    let (sum0, carry1) = sum0.overflowing_add(tmp2.wrapping_shl(64));
    let sum1 = tmp3
        .wrapping_add(tmp1.wrapping_shr(64))
        .wrapping_add(tmp2.wrapping_shr(64))
        .wrapping_add(carry0 as u128)
        .wrapping_add(carry1 as u128);

    (sum0, sum1)
}
/// Calculates c = a * b, without the possibility to overflow.
#[inline]
pub(crate) fn mul(a: &[u128; 2], b: &[u128; 2]) -> [u128; 4] {
    let mut res = [0; 2 * 2];
    let mut carry;

    (res[0], carry) = widening_mul(a[0], b[0]);
    (res[1], res[2]) = carrying_mul(a[0], b[1], carry);

    (res[1], carry) = carrying_mul(a[1], b[0], res[1]);
    (res[2], res[3]) = mac_with_carry(res[2], a[1], b[1], carry);

    res
}
/// Calculates the square lhs^2 without the possibility to overflow.
/// This returns the low-order (wrapping) bits and the high-order (overflow) bits of the result as two separate values, in that order.
/// This is faster than widening_mul(lhs, lhs)
#[inline]
pub(crate) fn widening_square(lhs: u128) -> (u128, u128) {
    let (lhs_hi, lhs_lo) = split_u128(lhs);
    let tmp0 = lhs_lo.wrapping_mul(lhs_lo);
    let tmp1 = lhs_lo.wrapping_mul(lhs_hi);
    let tmp3 = lhs_hi.wrapping_mul(lhs_hi);
    // tmp1 and tmp2 straddle the boundary. We have to handle three carries
    let (sum0, carry0) = tmp0.overflowing_add(tmp1.wrapping_shl(65));
    let sum1 = tmp3
        .wrapping_add(tmp1.wrapping_shr(63))
        .wrapping_add(carry0 as u128);

    (sum0, sum1)
}
/// Calculates the square lhs^2 + carry without the possibility to overflow.
/// This returns the low-order (wrapping) bits and the high-order (overflow) bits of the result as two separate values, in that order.
/// This is faster than carrying_mul(lhs, lhs, carry)
#[inline]
pub(crate) fn carrying_square(lhs: u128, carry: u128) -> (u128, u128) {
    let (sum0, sum1) = widening_square(lhs);
    let (sum0, carry2) = sum0.overflowing_add(carry);
    let sum1 = sum1.wrapping_add(carry2 as u128);

    (sum0, sum1)
}
/// Calculates c = a^2, without the possibility to overflow.
#[inline]
pub(crate) fn square(a: &[u128; 2]) -> [u128; 4] {
    let mut res = [0; 2 * 2];
    let mut carry;

    (res[0], res[1]) = widening_square(a[0]);
    (res[2], res[3]) = widening_square(a[1]);

    let (r1, r2) = widening_mul(a[0], a[1]);

    (res[1], carry) = res[1].overflowing_add(r1.wrapping_shl(1));
    (res[2], carry) = carrying_add(res[2], r2.wrapping_shl(1) | r1.wrapping_shr(127), carry);
    (res[3], _) = carrying_add(res[3], r2.wrapping_shr(127), carry);

    res
}
#[inline]
pub(crate) fn reduce(res: &mut [u128; 4], modulus: &[u128; 2], inv: u128) -> bool {
    let mut carry;

    let k = inv.wrapping_mul(res[0]);
    (_, carry) = carrying_mul(k, modulus[0], res[0]);
    (res[1], carry) = mac_with_carry(res[1], k, modulus[1], carry);
    res[2] = adc(res[2], 0, &mut carry);
    let r3 = carry;

    let k = inv.wrapping_mul(res[1]);
    (_, carry) = carrying_mul(k, modulus[0], res[1]);
    (res[2], carry) = mac_with_carry(res[2], k, modulus[1], carry);
    res[3] = adc(res[3], r3, &mut carry);

    carry != 0
}
/// Calculates a + b + carry and returns a tuple containing the sum and the output carry.
// u128::carrying_add is a nightly-only experimental API now. so use this version firstly.
// come from core/src/num/uint_macros.rs
#[inline]
pub(crate) fn carrying_add(a: u128, b: u128, carry: bool) -> (u128, bool) {
    let (r, c) = a.overflowing_add(b);
    let (x, d) = r.overflowing_add(carry as u128);
    (x, c || d)
}
/// Calculates a − b − borrow and returns a tuple containing the difference and the output borrow.
// u128::borrowing_sub is a nightly-only experimental API now. so use this version firstly.
#[inline]
pub(crate) fn borrowing_sub(a: u128, b: u128, borrow: bool) -> (u128, bool) {
    let (r, c) = a.overflowing_sub(b);
    let (x, d) = r.overflowing_sub(borrow as u128);
    (x, c || d)
}
//a + b * c + carry,
// This returns the low-order (wrapping) bits and the high-order (overflow) bits of the result as two separate values, in that order.
#[inline]
pub(crate) fn mac_with_carry(a: u128, b: u128, c: u128, carry: u128) -> (u128, u128) {
    let (lo, hi) = carrying_mul(b, c, a);
    let (r0, c) = lo.overflowing_add(carry);
    (r0, hi.wrapping_add(c as u128))
}
