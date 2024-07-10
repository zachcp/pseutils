use pymol_session_utils::{load_and_print_pickle, PSEData};
use serde_pickle::{from_reader, DeOptions};
use std::fs::File;
use std::io::Read;

#[test]
fn test_load_pickle_file() {
    let result = load_and_print_pickle();
    assert!(result.is_ok());
}

#[test]
fn test_deserialize_pickle_data() {
    let mut file = File::open("examples/example_molecule_only.pse").unwrap();
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).unwrap();

    let options = DeOptions::new()
        .replace_unresolved_globals()
        .decode_strings();

    let deserialized: PSEData = from_reader(&buffer[..], options).unwrap();

    // Add assertions to check the deserialized data
    // For example:
    // assert_eq!(deserialized.some_field, expected_value);
}

#[test]
fn test_pse_data_structure() {
    // Test the structure of your PSEData
    // For example:
    // let pse_data = PSEData { /* initialize with test data */ };
    // assert_eq!(pse_data.some_field, expected_value);
}

#[test]
fn test_specific_values() {
    let mut file = File::open("examples/example_molecule_only.pse").unwrap();
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).unwrap();

    let options = DeOptions::new()
        .replace_unresolved_globals()
        .decode_strings();

    let deserialized: PSEData = from_reader(&buffer[..], options).unwrap();

    // Add assertions to check specific values in the deserialized data
    // For example:
    // assert_eq!(deserialized.main[0].some_field, expected_value);
}

#[test]
fn test_error_handling() {
    // Test error handling for various scenarios
    // For example:
    // - Test with a non-existent file
    // - Test with an invalid pickle file
    // - Test with unexpected data structure
}

#[test]
fn test_performance() {
    // Add performance tests if needed
    // For example, measure the time it takes to load and deserialize a large file
}
