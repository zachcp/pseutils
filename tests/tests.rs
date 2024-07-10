use pymol_session_utils::PSEData;

#[test]
fn test_load_pse_data() {
    let deserialized: PSEData = PSEData::load("tests/data/example_molecule_only.pse").unwrap();
    // Check fields of PSEData
    assert!(deserialized.version == 3000000);
}
