# rs-asat-hor
Library for parsing and iterating through alpha-satellite higher-order repeats.

### Usage
```rust
use rs_asat_hor::{HOR, Monomer};

let mons = [
    Monomer::new("S1C1/5/19H1L.1").unwrap(),
    Monomer::new("S1C1/5/19H1L.2").unwrap(),
    Monomer::new("S1C1/5/19H1L.3").unwrap(),
];
let hor = HOR::from_monomers(&mons).unwrap();
assert_eq!(
    format!("{}", hor[0]),
    "S1C1/5/19H1L.1-3"
)
```

To add this crate:
```bash
cargo add --git https://github.com/koisland/rs-asat-hor/tree/main rs-asat-hor
```

### Why?
* All existing HOR stv code is string-y.
* Edge cases.
* Nothing able to go from stv to monomers.

### TODO
* Python binding.
* Documentation
