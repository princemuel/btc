#![allow(unused)]

pub const fn mulmod(mut a: u128, mut b: u128, m: u128) -> u128 {
    a = a.rem_euclid(m);
    b = b.rem_euclid(m);

    let mut result = 0u128;

    while b > 0 {
        if b & 1 == 1 {
            result = addmod(result, a, m);
        }
        a = addmod(a, a, m);
        b >>= 1;
    }

    result
}

pub const fn addmod(mut a: u128, mut b: u128, m: u128) -> u128 {
    let s = a + b;

    if s >= m { s - m } else { s }

    // match a.checked_add(b) {
    //     Some(s) => {
    //         if s >= m {
    //             s - m
    //         } else {
    //             s
    //         }
    //     },
    //     None => {
    //         // Handle overflow in addition
    //         let a_mod = a.rem_euclid(m);
    //         let b_mod = b.rem_euclid(m);
    //         a_mod.wrapping_add(b_mod).rem_euclid(m)
    //     },
    // }
}

pub const fn powmod(mut a: i128, mut b: u128, m: i128) -> i128 {
    let mut result = 1i128;

    while b > 0 {
        if b & 1 == 1 {
            result = (result * a).rem_euclid(m);
        }
        a = (a * a).rem_euclid(m);
        b >>= 1;
    }
    result
}

// pub const fn invmod()

/// Extended Euclidean Algorithm
/// Returns (gcd, x, y) where ax + by = gcd(a, b)
pub fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if a == 0 {
        let (gcd, y) = if b >= 0 { (b, 1) } else { (-b, -1) };
        return (gcd, 0, y);
    }

    let r = b.rem_euclid(a);
    let q = b.div_euclid(a);

    let (gcd, x1, y1) = extended_gcd(r, a);

    let x = y1 - q * x1;
    let y = x1;

    (gcd, x, y)
}

pub fn step(a: &mut i128, old_a: &mut i128, quotient: i128) {
    let temp = *a;
    *a = *old_a - quotient * temp;
    *old_a = temp;
}

pub fn extended_euclidean_algorithm(a: i128, b: i128) -> (i128, i128, i128) {
    // Remainders
    let (mut old_r, mut rem) = (a, b);
    // Bézout coefficients for a
    let (mut old_s, mut coeff_s) = (1, 0);
    // Bézout coefficients for b
    let (mut old_t, mut coeff_t) = (0, 1);

    while rem != 0 {
        let quotient = old_r.div_euclid(rem);

        // r_k+1 = r_k-1 - q_k * r_k
        step(&mut rem, &mut old_r, quotient);
        // s_k+1 = s_k-1 - q_k * s_k
        step(&mut coeff_s, &mut old_s, quotient);
        // t_k+1 = t_k-1 - q_k * t_k
        step(&mut coeff_t, &mut old_t, quotient);
    }

    if old_r < 0 {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
    }

    (old_r, old_s, old_t)
}

#[inline]
pub fn gcd(a: i128, b: i128) -> i128 {
    let (g, ..) = extended_euclidean_algorithm(a, b);
    g
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};

    use super::*;

    fn mod_inverse(a: i128, m: i128) -> i128 {
        let (g, x, _) = extended_euclidean_algorithm(a, m);
        assert_eq!(g, 1, "inverse does not exist");
        x.rem_euclid(m)
    }

    // Helper: brute-force mod inverse for cross-checking (only for small m)
    fn mod_inverse_bruteforce(a: i128, m: i128) -> Option<i128> {
        (0..m).find(|&x| (a * x).rem_euclid(m) == 1)
    }

    #[test]
    fn test_extended_gcd_basic() {
        assert_eq!(extended_euclidean_algorithm(0, 7), (7, 0, 1));
        assert_eq!(extended_euclidean_algorithm(7, 0), (7, 1, 0));

        assert_eq!(extended_euclidean_algorithm(7, 13).0, 1);
        assert_eq!(extended_euclidean_algorithm(13, 7).0, 1);

        assert_eq!(extended_euclidean_algorithm(48, 18).0, 6);
        assert_eq!(extended_euclidean_algorithm(18, 48).0, 6);

        assert_eq!(extended_euclidean_algorithm(270, 192).0, 6);
        assert_eq!(extended_euclidean_algorithm(192, 270).0, 6);
    }

    #[test]
    fn test_extended_gcd_bezout_identity() {
        let range = -50..=50;

        for a in range.clone() {
            for b in range.clone() {
                let (g, x, y) = extended_euclidean_algorithm(a, b);
                // gcd must match reference implementation
                assert_eq!(g, gcd(a, b), "gcd mismatch for a={a}, b={b}");

                // gcd must be non negative
                assert!(g >= 0, "gcd must be non-negative");

                // Bezout identity must hold
                assert_eq!(x * a + y * b, g, "Bezout identity failed for a={a}, b={b}");
            }
        }
    }

    #[test]
    fn test_extended_gcd_randomized() {
        use rand::Rng;
        let mut rng = rand::rng();

        for _ in 0..10_000 {
            let a = rng.random_range(-10_000..10_000);
            let b = rng.random_range(-10_000..10_000);

            let (g, x, y) = extended_euclidean_algorithm(a, b);

            // Check Bézout identity
            assert_eq!(x * a + y * b, g);

            // Check correctness of gcd sign
            if a == 0 && b == 0 {
                assert_eq!(g, 0);
            } else {
                assert!(g >= 0);
            }
        }
    }

    #[test]
    fn test_mod_inverse_small_modulus() {
        let m = 7;

        let table = |m: i128| {
            assert!(m > 1, "modulus must be greater than 1");
            (0..m)
                .filter_map(|a| {
                    let (g, x, _) = extended_euclidean_algorithm(a, m);
                    (g == 1).then(|| (a, x.rem_euclid(m)))
                })
                .collect::<Vec<_>>()
        };

        // Convert to a map for convenience
        let map: HashMap<_, _> = table(m).iter().copied().collect();

        // 1. Check correctness: a * inv ≡ 1 (mod m)
        for (&a, &inv) in &map {
            assert_eq!(
                (a * inv).rem_euclid(m),
                1,
                "incorrect inverse: {a} * {inv} mod {m} != 1",
            );
        }

        // 2. Check that inverse is involutive: inv(inv(a)) = a
        for (&a, &inv) in &map {
            let inv2 = map[&inv];
            assert_eq!(inv2, a, "involution failed: inv(inv({a})) = {inv2} != {a}");
        }

        // 3. Check completeness: number of invertible residues = φ(m)
        let phi = (1..m).filter(|&x| gcd(x, m) == 1).count();
        let map_len = map.len();
        assert_eq!(map_len, phi, "inverse table size {map_len} != φ({m}) = {phi}",);

        // 4. Check uniqueness: every inverse is unique
        let values: HashSet<_> = map.values().copied().collect();
        assert_eq!(values.len(), map.len(), "modular inverses are not unique");
    }

    #[test]
    fn test_mod_inverse_matches_bruteforce() {
        for m in 2..50 {
            for a in -200..200 {
                let fast = extended_euclidean_algorithm(a, m);
                if fast.0 == 1 {
                    let inv1 = fast.1.rem_euclid(m);
                    let inv2 = mod_inverse_bruteforce(a, m).unwrap();
                    assert_eq!(inv1, inv2, "m={}, a={}", m, a);
                }
            }
        }
    }

    #[test]
    fn test_mod_inverse_large_primes() {
        // Prime moduli
        let primes = [1_000_000_007i128, 1_000_000_009i128, 2_147_483_647i128];

        for &p in &primes {
            for a in [1, 2, 3, 12345, -7, -999_999] {
                let inv = mod_inverse(a, p);
                assert_eq!((inv * a).rem_euclid(p), 1);
            }
        }
    }

    #[test]
    fn test_mod_inverse_random() {
        use rand::Rng;
        let mut rng = rand::rng();

        for _ in 0..10_000 {
            let m = rng.random_range(2..100_000);
            let a = rng.random_range(-100_000..100_000);

            let (g, x, _) = extended_euclidean_algorithm(a, m);

            if g == 1 {
                let inv = x.rem_euclid(m);
                assert_eq!((inv * a).rem_euclid(m), 1);
            }
        }
    }
}
