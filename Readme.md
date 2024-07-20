# pymol_session_utils

**Note: Work In Progress. Feature Incomplete**


## Goal

Binary serializaton from `.pse` to `X` using RUST. Where I want `X` to be [mol-view-spec](https://github.com/molstar/mol-view-spec)

## Testing

```sh
cargo test

# to get all the print statements that are sprinkled around.
cargo test -- --nocapture

# build the docs
cargo clean --doc && cargo doc --no-deps --open

# to look at results
cargo test # generate an example
cd test_temporary && python -m http.server
```

## Status

- PSE Conversion to Molstar-ready formula.
- Basic processing of molecules and selections

![](resources/images/pymol_molstar.png)

## Related

- [pymol source](https://github.com/schrodinger/pymol-open-source)
- [mol-view-spec](https://github.com/molstar/mol-view-spec)
- [MichelaNGLo-transpiler](https://github.com/matteoferla/MichelaNGLo-transpiler)
  - [conversion notes](https://github.com/matteoferla/MichelaNGLo-transpiler/blob/master/docs/conversion.md)
  - [view rotation notes](https://github.com/matteoferla/MichelaNGLo-transpiler/blob/master/docs/notes_on_view_conversion.md)
