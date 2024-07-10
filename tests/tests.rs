use pymol_session_utils::PSEData;

#[test]
fn test_load_pse_data_molecule_only() {
    let deserialized: PSEData = PSEData::load("tests/data/example_molecule_only.pse").unwrap();
    // deserialized.to_json("tests/data/example_molecule_only.json");
    // Check fields of PSEData
    assert!(deserialized.version == 3000000);
}
