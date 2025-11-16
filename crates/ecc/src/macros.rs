macro_rules! impl_from_int_for_primefield_signed {
    ($($t:ty),*) => {
        $(
            impl<const P: u128> From<$t> for PrimeField<P> {
                #[inline(always)]
                fn from(value: $t) -> Self {
                    // Signed: rem_euclid handles negatives properly
                    Self((value as i128).rem_euclid(P as i128))
                }
            }
        )*
    };
}
pub(crate) use impl_from_int_for_primefield_signed;

macro_rules! impl_from_int_for_primefield_unsigned {
    ($($t:ty),*) => {
        $(
            impl<const P: u128> From<$t> for PrimeField<P> {
                #[inline(always)]
                fn from(value: $t) -> Self {
                    Self(((value as u128) % P) as i128)
                }
            }
        )*
    };
}
pub(crate) use impl_from_int_for_primefield_unsigned;

macro_rules! impl_tryfrom_primefield_for_small_ints {
    ($($t:ty),*) => {
        $(
            impl<const P: u128> TryFrom<PrimeField<P>> for $t {
                type Error = &'static str;

                #[inline(always)]
                fn try_from(value: PrimeField<P>) -> Result<Self, Self::Error> {
                    <$t>::try_from(value.0)
                        .map_err(|_| "PrimeField element does not fit in target type")
                    }
            }
        )*
    };
}

pub(crate) use impl_tryfrom_primefield_for_small_ints;
