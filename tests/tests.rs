use pymol_session_utils::psedata::PSEData;

#[test]
fn test_load_pse_data_molecule_only() {
    //https://users.rust-lang.org/t/serde-untagged-enum-ruins-precise-errors/54128/2
    // https://www.gustavwengel.dk/serde-untagged-enum-errors-are-bad
    let deserialized: PSEData = PSEData::load("tests/data/example_molecule_only.pse").unwrap();
    // deserialized.to_json("tests/data/example_molecule_only.json");
    // Check fields of PSEData
    assert!(deserialized.version == 3000000);
}

#[test]
fn test_load_pse_data_molecule_selection() {
    let deserialized: PSEData = PSEData::load("tests/data/example.pse").unwrap();
    assert!(deserialized.version == 3000000);
}
