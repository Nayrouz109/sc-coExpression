

def list_keys_and_subkeys(d, parent_key=""):
    """
    Recursively lists all keys and sub-keys of a nested dictionary.

    Args:
        d (dict): The dictionary to process.
        parent_key (str): Used for nested keys; defaults to an empty string.

    Returns:
        None: Prints keys and sub-keys.
    """
    for key, value in d.items():
        # Combine parent key with the current key for better visualization.
        full_key = f"{parent_key}.{key}" if parent_key else key
        print(full_key)  # Print the full key path
        
        # Check if the value is a nested dictionary and recurse.
        if isinstance(value, dict):
            list_keys_and_subkeys(value, full_key)


# print out the length of dictionary components 
for key, value in resultsv1["pairwise_intersections"].items():
    print(f"Key: {key}, Length: {len(value)}")


# Example dictionary
nested_dict = {
    "cell_type1": {
        "study1": ["geneA", "geneB"],
        "study2": ["geneC"],
    },
    "cell_type2": {
        "study3": ["geneD", "geneE"],
        "study4": ["geneF"],
    }
}

# Call the function
list_keys_and_subkeys(nested_dict)
