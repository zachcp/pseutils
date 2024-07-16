# pymol_session_utils

**WIP: Currently not working**

## Testing

```sh
cargo test

# to get all the print statements that are sprinkled around.
cargo test -- --nocapture
```


## Goal

Binary serializaton from `.pse` to `X` using RUST. Where I want `X` to be [mol-view-spec](https://github.com/molstar/mol-view-spec)

## Related

- [pymol source](https://github.com/schrodinger/pymol-open-source)
- [mol-view-spec](https://github.com/molstar/mol-view-spec)
- [MichelaNGLo-transpiler](https://github.com/matteoferla/MichelaNGLo-transpiler)
  - [conversion notes](https://github.com/matteoferla/MichelaNGLo-transpiler/blob/master/docs/conversion.md)
  - [view rotation notes](https://github.com/matteoferla/MichelaNGLo-transpiler/blob/master/docs/notes_on_view_conversion.md)
