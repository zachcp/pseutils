use pymol_session_utils::PSEData;

#[test]
fn test_load_pse_data_molecule_only() {
    let deserialized: PSEData = PSEData::load("tests/data/example_molecule_only.pse").unwrap();
    // Check fields of PSEData
    assert!(deserialized.version == 3000000);
}

// #[test]
// fn test_load_pse_data_molecule_selection() {
//     let deserialized: PSEData = PSEData::load("tests/data/example.pse").unwrap();
//     // Check fields of PSEData
//     assert!(deserialized.version == 3000000);
// }
