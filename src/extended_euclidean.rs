#[derive(Debug, PartialEq)]
pub struct ExtendedEuclideanResult {
    pub gcd: i64,
    pub bezout_coefficient_a: i64,
    pub bezout_coefficient_b: i64,
}

impl ExtendedEuclideanResult {
    pub fn new(gcd: i64, bezout_coefficient_a: i64, bezout_coefficient_b: i64) -> Self {
        Self {
            gcd,
            bezout_coefficient_a,
            bezout_coefficient_b,
        }
    }
}

/// Runs the extended euclidean algorithm on `a` and `b`.
pub fn extended_euclidean(a: i64, b: i64) -> ExtendedEuclideanResult {
    let mut remainder_prev;
    let mut remainder;
    let mut s_prev;
    let mut s;
    let mut t_prev;
    let mut t;

    if a > b {
        remainder_prev = a;
        remainder = b;
        s_prev = 1;
        s = 0;
        t_prev = 0;
        t = 1;
    } else {
        remainder_prev = b;
        remainder = a;
        s_prev = 0;
        s = 1;
        t_prev = 1;
        t = 0;
    }

    while remainder != 0 {
        let quotient = remainder_prev / remainder;

        let previous_remainder = remainder;
        remainder = remainder_prev - quotient * remainder;
        remainder_prev = previous_remainder;

        let previous_s = s;
        s = s_prev - quotient * s;
        s_prev = previous_s;

        let previous_t = t;
        t = t_prev - quotient * t;
        t_prev = previous_t;
    }

    // Return GCD as a nonnegative integer, adjusting the coefficients accordingly.
    if remainder_prev < 0 {
        remainder_prev = remainder_prev * -1;
        s_prev = s_prev * -1;
        t_prev = t_prev * -1;
    }

    ExtendedEuclideanResult::new(remainder_prev, s_prev, t_prev)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eea() {
        let tests = [
            (240, 46, 2),
            (46, 240, 2),
            (-240, 46, 2),
            (-46, 240, 2),
            (46, -240, 2),
            (0, 0, 0),
            (1, 0, 1),
            (0, 1, 1),
            (1030203, 4393920, 3),
            (92874881, 2343483, 1),
        ];

        for (a, b, expected_gcd) in tests {
            let result = extended_euclidean(a, b);
            assert_eq!(result.gcd, expected_gcd);
            assert_eq!(
                a * result.bezout_coefficient_a + b * result.bezout_coefficient_b,
                expected_gcd,
            );
        }
    }
}
