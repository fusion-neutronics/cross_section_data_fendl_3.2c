#!/usr/bin/env python3
# Script to download and extract FENDL 3.2c nuclear data
# before runnung this script
# pip install openmc_data
# install njoy and openmc
# generate_fendl -r 3.2c

import glob
import json
import sys
import tarfile
from pathlib import Path
from typing import Optional

import numpy as np
import openmc

# Import monkey patches to add simplified methods to OpenMC classes
from openmc_patches import *


class CompactJSONEncoder(json.JSONEncoder):
    """A JSON Encoder that puts arrays of primitives on a single line."""

    def __init__(self, *args, **kwargs):
        # Default to 2-space indentation if not specified
        if kwargs.get("indent") is None:
            kwargs["indent"] = 2
        # Ensure compact separators (no spaces) for arrays
        if kwargs.get("separators") is None:
            kwargs["separators"] = (',', ':')
        super().__init__(*args, **kwargs)
        self.indentation_level = 0

    def default(self, o):
        """Handle NumPy types"""
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        return json.JSONEncoder.default(self, o)

    def encode(self, o):
        """Encode JSON objects with special handling for arrays of primitives."""
        if isinstance(o, (list, tuple)):
            # Put arrays of primitives on a single line
            if all(isinstance(x, (int, float, bool, str, type(None))) for x in o):
                return "[" + ",".join(json.dumps(el) for el in o) + "]"
            # For other containers, use standard formatting
            return self._encode_container(o, "[", "]")
        if isinstance(o, dict):
            if not o:
                return "{}"
            # Ensure keys are strings
            o = {str(k) if k is not None else "null": v for k, v in o.items()}
            if self.sort_keys:
                o = dict(sorted(o.items(), key=lambda x: x[0]))
            # Format dictionary
            return self._encode_container(o, "{", "}")
        # For primitive types, use standard JSON encoding
        return json.dumps(
            o,
            skipkeys=self.skipkeys,
            ensure_ascii=self.ensure_ascii,
            check_circular=self.check_circular,
            allow_nan=self.allow_nan,
            sort_keys=self.sort_keys,
            separators=(self.item_separator, self.key_separator),
            default=self.default,
        )

    def _encode_container(self, container, open_char, close_char):
        """Encode a container with proper indentation."""
        is_dict = isinstance(container, dict)

        # Handle very small containers on single line (dicts with 1 item or less)
        if is_dict and len(container) <= 1:
            if not container:
                return "{}"
            key, value = next(iter(container.items()))
            return f"{{ {self.encode(key)}: {self.encode(value)} }}"

        # Format container across multiple lines
        self.indentation_level += 1
        if is_dict:
            items = [
                f"{self.indent_str}{self.encode(k)}: {self.encode(v)}"
                for k, v in container.items()
            ]
        else:
            items = [f"{self.indent_str}{self.encode(el)}" for el in container]
        self.indentation_level -= 1

        return f"{open_char}\n" + ",\n".join(items) + f"\n{self.indent_str}{close_char}"

    def iterencode(self, data, **kwargs):
        """Required to also work with `json.dump`."""
        return self.encode(data)

    @property
    def indent_str(self) -> str:
        """Get the indentation string for the current level."""
        if isinstance(self.indent, int):
            return " " * (self.indentation_level * self.indent)
        elif isinstance(self.indent, str):
            return self.indentation_level * self.indent
        else:
            return ""


def extract_cross_section_data(
    h5_filename: str,
    library: str,
    json_filename: Optional[str] = None,
) -> dict:
    """Extract cross section data and scattering data from HDF5 and write to JSON using OpenMC monkey patches."""
    
    # Load the isotope data
    isotope_object = openmc.data.IncidentNeutron.from_hdf5(h5_filename)
    isotope_object.library = library
    
    # Use monkey patch to extract cross-section data with hierarchical MTs
    result = isotope_object.add_hierarchical_mts()
    
    # Extract scattering data and integrate into each reaction
    print(f"Extracting scattering data for {isotope_object.name}...")
    scattering_count = 0
    
    for reaction_mt, reaction in isotope_object.reactions.items():
        if reaction.products:
            neutron_products = []
            
            for product in reaction.products:
                if product.particle == 'neutron':
                    try:
                        # Use monkey patch to convert product to dictionary
                        product_dict = product.to_dict()
                        neutron_products.append(product_dict)
                    except Exception as e:
                        print(f"  Warning: Could not extract product data for MT {reaction_mt}: {e}")
            
            if neutron_products:
                # Add scattering data to each temperature for this reaction
                for temperature in result["reactions"]:
                    if str(reaction_mt) in result["reactions"][temperature]:
                        result["reactions"][temperature][str(reaction_mt)]["products"] = neutron_products
                scattering_count += 1
    
    if scattering_count > 0:
        print(f"  Found scattering data for {scattering_count} reactions")
    else:
        print("  No scattering data found")
    
    # Set output filename
    if json_filename is None:
        json_filename = f"{isotope_object.name}.json"
    
    # Sort MTs numerically in each temperature before writing
    for temperature in result["reactions"]:
        reactions = result["reactions"][temperature]
        sorted_reactions = dict(sorted(reactions.items(), key=lambda x: int(x[0])))
        result["reactions"][temperature] = sorted_reactions
    
    # Write to JSON file using custom encoder
    with open(json_filename, "w") as f:
        json.dump(result, f, cls=CompactJSONEncoder, indent=2)
    
    return result


h5_files = glob.glob("fendl-3.2c-hdf5/*.h5")

for i, file in enumerate(h5_files):
        print(f"Processing file {i+1}/{len(h5_files)}: {file}")
        extract_cross_section_data(h5_filename=file, library="FENDL-3.2C")
