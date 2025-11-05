"""
Monkey patches for OpenMC classes to add JSON/dict serialization capabilities.
Import this module to automatically add to_dict() methods to OpenMC classes.
"""

import openmc.data
import json
import numpy as np

# Create reverse mapping from element symbol to element name for efficiency
SYMBOL_TO_ELEMENT = {symbol: name for name, symbol in openmc.data.ELEMENT_SYMBOL.items()}

def _convert_to_list(data):
    """Helper function to convert numpy arrays or iterables to lists"""
    if data is None:
        return None
    if hasattr(data, 'tolist'):
        return data.tolist()
    elif hasattr(data, '__iter__') and not isinstance(data, str):
        try:
            # Handle nested structures recursively
            return [_convert_to_list(item) for item in data]
        except:
            return list(data)
    else:
        return data

def _convert_distribution_list(dist_list):
    """Helper function to convert a list of distributions to dictionary format"""
    if dist_list is None:
        return None
    
    result = []
    for dist in dist_list:
        if hasattr(dist, 'to_dict'):
            # If the distribution has its own to_dict method
            result.append(dist.to_dict())
        elif hasattr(dist, 'x') and hasattr(dist, 'p'):
            # Tabulated distribution
            result.append({
                'type': 'Tabulated',
                'x': _convert_to_list(dist.x),
                'p': _convert_to_list(dist.p)
            })
        elif hasattr(dist, 'parameters'):
            # Parametric distribution
            result.append({
                'type': type(dist).__name__,
                'parameters': dist.parameters
            })
        else:
            result.append({
                'type': type(dist).__name__,
                'raw': str(dist)
            })
    return result

# Monkey patch UncorrelatedAngleEnergy
def uncorrelated_angle_energy_to_dict(self):
    """Convert UncorrelatedAngleEnergy to dictionary representation"""
    result = {
        'type': 'UncorrelatedAngleEnergy'
    }
    
    # Handle angle distribution
    if hasattr(self, 'angle') and self.angle is not None:
        result['angle'] = self.angle.to_dict()
        
    # Handle energy distribution  
    if hasattr(self, 'energy') and self.energy is not None:
        result['energy'] = self.energy.to_dict()
        
    return result

openmc.data.UncorrelatedAngleEnergy.to_dict = uncorrelated_angle_energy_to_dict

# Monkey patch CorrelatedAngleEnergy
def correlated_angle_energy_to_dict(self):
    """Convert CorrelatedAngleEnergy to dictionary representation"""
    result = {
        'type': 'CorrelatedAngleEnergy'
    }
    
    if hasattr(self, 'energy') and self.energy is not None:
        result['energy'] = _convert_to_list(self.energy)
        
    if hasattr(self, 'energy_out') and self.energy_out is not None:
        result['energy_out'] = _convert_distribution_list(self.energy_out)
        
    if hasattr(self, 'mu') and self.mu is not None:
        result['mu'] = []
        for mu_group in self.mu:
            if hasattr(mu_group, '__iter__'):
                result['mu'].append(_convert_distribution_list(mu_group))
            else:
                result['mu'].append(str(mu_group))
                
    if hasattr(self, 'breakpoints') and self.breakpoints is not None:
        result['breakpoints'] = self.breakpoints
    if hasattr(self, 'interpolation') and self.interpolation is not None:
        result['interpolation'] = _convert_to_list(self.interpolation)
        
    return result

openmc.data.CorrelatedAngleEnergy.to_dict = correlated_angle_energy_to_dict

# Monkey patch AngleDistribution (base class for angle distributions)
def angle_distribution_to_dict(self):
    """Convert angle distribution to dictionary representation"""
    result = {
        'type': type(self).__name__
    }
    
    if hasattr(self, 'energy') and self.energy is not None:
        result['energy'] = _convert_to_list(self.energy)
        
    if hasattr(self, 'mu') and self.mu is not None:
        result['mu'] = _convert_distribution_list(self.mu)
        
    if hasattr(self, 'breakpoints') and self.breakpoints is not None:
        result['breakpoints'] = self.breakpoints
    if hasattr(self, 'interpolation') and self.interpolation is not None:
        result['interpolation'] = _convert_to_list(self.interpolation)
        
    return result

openmc.data.AngleDistribution.to_dict = angle_distribution_to_dict

# Monkey patch EnergyDistribution (base class for energy distributions) 
def energy_distribution_to_dict(self):
    """Convert energy distribution to dictionary representation"""
    result = {
        'type': type(self).__name__
    }
    
    if hasattr(self, 'energy') and self.energy is not None:
        result['energy'] = _convert_to_list(self.energy)
        
    if hasattr(self, 'energy_out') and self.energy_out is not None:
        result['energy_out'] = _convert_distribution_list(self.energy_out)
        
    if hasattr(self, 'parameters') and self.parameters is not None:
        result['parameters'] = self.parameters
        
    if hasattr(self, 'breakpoints') and self.breakpoints is not None:
        result['breakpoints'] = self.breakpoints
    if hasattr(self, 'interpolation') and self.interpolation is not None:
        result['interpolation'] = _convert_to_list(self.interpolation)
        
    return result

openmc.data.EnergyDistribution.to_dict = energy_distribution_to_dict

# Monkey patch KalbachMann
def kalbach_mann_to_dict(self):
    """Convert KalbachMann to dictionary representation"""
    result = {
        'type': 'KalbachMann'
    }
    
    if hasattr(self, 'energy') and self.energy is not None:
        result['energy'] = _convert_to_list(self.energy)
        
    if hasattr(self, 'energy_out') and self.energy_out is not None:
        result['energy_out'] = _convert_distribution_list(self.energy_out)
        
    if hasattr(self, 'precompound_fraction') and self.precompound_fraction is not None:
        result['precompound_fraction'] = _convert_distribution_list(self.precompound_fraction)
        
    if hasattr(self, 'slope') and self.slope is not None:
        result['slope'] = _convert_distribution_list(self.slope)
        
    return result

openmc.data.KalbachMann.to_dict = kalbach_mann_to_dict

# Monkey patch Tabulated1D
def tabulated1d_to_dict(self):
    """Convert Tabulated1D to dictionary representation"""
    result = {
        'type': 'Tabulated1D',
        'x': _convert_to_list(self.x),
        'y': _convert_to_list(self.y),
        'breakpoints': self.breakpoints,
        'interpolation': _convert_to_list(self.interpolation)
    }
    return result

openmc.data.Tabulated1D.to_dict = tabulated1d_to_dict

# Monkey patch Polynomial
def polynomial_to_dict(self):
    """Convert Polynomial to dictionary representation"""
    result = {
        'type': 'Polynomial',
        'coefficients': _convert_to_list(self.coef)
    }
    return result

openmc.data.Polynomial.to_dict = polynomial_to_dict

# Monkey patch Product
def product_to_dict(self):
    """Convert Product to dictionary representation"""
    result = {
        'particle': self.particle,
        'emission_mode': self.emission_mode,
        'decay_rate': self.decay_rate,
        'applicability': [],
        'distribution': []
    }
    
    # Handle applicability
    for app in self.applicability:
        if hasattr(app, 'to_dict'):
            result['applicability'].append(app.to_dict())
        else:
            result['applicability'].append({
                'type': type(app).__name__,
                'raw': str(app)
            })
    
    # Handle distribution
    for dist in self.distribution:
        if hasattr(dist, 'to_dict'):
            result['distribution'].append(dist.to_dict())
        else:
            result['distribution'].append({
                'type': type(dist).__name__,
                'raw': str(dist)
            })
    
    # Handle yield - ensure it's JSON serializable
    if hasattr(self.yield_, 'to_dict'):
        result['yield'] = self.yield_.to_dict()
    else:
        result['yield'] = {
            'type': type(self.yield_).__name__,
            'values': _convert_to_list(self.yield_)
        }
    
    return result

openmc.data.Product.to_dict = product_to_dict

# Add convenience method to convert products list to JSON
def products_to_json(products, indent=2):
    """Convert a list of products to JSON string"""
    products_dict = [product.to_dict() for product in products]
    return json.dumps(products_dict, indent=indent)

# Custom JSON encoder to handle numpy arrays with compact formatting
class NumpyEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        # Use compact separators (no spaces) by default
        if 'separators' not in kwargs:
            kwargs['separators'] = (',', ':')
        super().__init__(*args, **kwargs)
    
    def default(self, obj):
        if hasattr(obj, 'tolist'):
            return obj.tolist()
        elif hasattr(obj, '__array__'):
            return obj.tolist()
        return super().default(obj)

# Make the function available as a module function
def to_json(obj, indent=2):
    """Convert any object with to_dict method to JSON string"""
    if hasattr(obj, 'to_dict'):
        return json.dumps(obj.to_dict(), indent=indent, cls=NumpyEncoder)
    elif hasattr(obj, '__iter__') and not isinstance(obj, str):
        # Handle lists of objects
        data = [item.to_dict() if hasattr(item, 'to_dict') else str(item) for item in obj]
        return json.dumps(data, indent=indent, cls=NumpyEncoder)
    else:
        return json.dumps(str(obj), indent=indent, cls=NumpyEncoder)

def safe_json_dump(data, file_path, indent=2):
    """Safely dump data to JSON file, handling numpy arrays with compact formatting"""
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=indent, cls=NumpyEncoder, separators=(',', ':'))

# Monkey patch IncidentNeutron for cross-section data extraction
def incident_neutron_to_cross_section_dict(self):
    """Convert IncidentNeutron to cross-section dictionary format for JSON export"""
    
    def normalize_temp_key(temp_key):
        """Convert temperature key to string without 'K' suffix"""
        return str(temp_key).replace("K", "")
    
    # Get element full name from symbol using OpenMC's ELEMENT_SYMBOL mapping
    element_name = SYMBOL_TO_ELEMENT.get(self.atomic_symbol, self.atomic_symbol.lower())
    
    result = {
        "atomic_number": self.atomic_number,
        "atomic_symbol": self.atomic_symbol,
        "atomic_weight_ratio": self.atomic_weight_ratio,
        "mass_number": self.mass_number,
        "metastable": self.metastable,
        "name": self.name,
        "element": element_name,
        "neutron_number": self.mass_number - self.atomic_number,
        "temperatures": [normalize_temp_key(temp) for temp in self.temperatures],
        "library": self.library,
    }
    
    # Collect all unique temperatures
    all_temperatures = set()
    for reaction in self.reactions.values():
        all_temperatures.update(normalize_temp_key(temp) for temp in reaction.xs.keys())
    
    # Initialize nested structure
    result["reactions"] = {}
    for temperature in all_temperatures:
        result["reactions"][temperature] = {}
    
    # Energy grid
    result["energy"] = {}
    for temperature, energy in self.energy.items():
        result["energy"][normalize_temp_key(temperature)] = energy
    
    # Extract reaction data
    for reaction_mt, reaction in self.reactions.items():
        for temperature, tabular in reaction.xs.items():
            temp_key = normalize_temp_key(temperature)
            result["reactions"][temp_key][str(reaction_mt)] = {
                "interpolation": tabular.interpolation,
                "threshold_idx": tabular._threshold_idx,
                "xs": tabular.y,
                "q_value": reaction.q_value
            }
    
    return result

openmc.data.IncidentNeutron.to_cross_section_dict = incident_neutron_to_cross_section_dict

# Monkey patch for direct JSON export
def incident_neutron_export_cross_sections_json(self, filename=None):
    """Export cross-section data directly to JSON file"""
    if filename is None:
        filename = f"{self.name}.json"
    
    data = self.to_cross_section_dict()
    
    # Sort MTs numerically in each temperature
    for temperature in data["reactions"]:
        reactions = data["reactions"][temperature]
        data["reactions"][temperature] = dict(sorted(reactions.items(), key=lambda x: int(x[0])))
    
    safe_json_dump(data, filename)
    return data

openmc.data.IncidentNeutron.export_cross_sections_json = incident_neutron_export_cross_sections_json

# Monkey patch to add MT sum rules capability
def incident_neutron_add_hierarchical_mts(self, sum_rules=None):
    """Add hierarchical MTs based on sum rules to the cross-section dictionary"""
    
    if sum_rules is None:
        sum_rules = {
            1: [2, 3],
            3: [4, 5, 11, 16, 17, 22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 35,
                36, 37, 41, 42, 44, 45, 152, 153, 154, 156, 157, 158, 159, 160,
                161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
                173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185,
                186, 187, 188, 189, 190, 194, 195, 196, 198, 199, 200],
            4: list(range(50, 92)),
            16: list(range(875, 892)),
            18: [19, 20, 21, 38],
            27: [18, 101],
            101: [102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114,
                  115, 116, 117, 155, 182, 191, 192, 193, 197],
            103: list(range(600, 650)),
            104: list(range(650, 700)),
            105: list(range(700, 750)),
            106: list(range(750, 800)),
            107: list(range(800, 850)),
        }
    
    data = self.to_cross_section_dict()
    
    # Apply sum rules to each temperature
    for temperature in data["reactions"]:
        reactions_at_temp = data["reactions"][temperature]
        available_mts = set(str(mt) for mt in reactions_at_temp.keys())
        energy_grid = data["energy"][temperature]
        grid_len = len(energy_grid)
        
        constructed = set()
        while True:
            added = False
            for mt, children in reversed(list(sum_rules.items())):
                mt_str = str(mt)
                if mt_str in reactions_at_temp:
                    continue
                    
                present_children = [str(child) for child in children if str(child) in available_mts]
                if not present_children:
                    continue
                
                # Construct the hierarchical MT
                aligned_arrays = []
                min_threshold_idx = grid_len
                
                for child in present_children:
                    child_xs = np.array(reactions_at_temp[child]["xs"])
                    child_threshold_idx = reactions_at_temp[child]["threshold_idx"]
                    min_threshold_idx = min(min_threshold_idx, child_threshold_idx)
                    
                    arr = np.zeros(grid_len)
                    arr[child_threshold_idx:child_threshold_idx + len(child_xs)] = child_xs
                    aligned_arrays.append(arr)
                
                summed_xs_full = np.sum(aligned_arrays, axis=0)
                summed_xs = summed_xs_full[min_threshold_idx:]
                interp = reactions_at_temp[present_children[0]]["interpolation"]
                
                reactions_at_temp[mt_str] = {
                    "interpolation": interp,
                    "threshold_idx": min_threshold_idx,
                    "xs": summed_xs.tolist(),
                }
                
                available_mts.add(mt_str)
                constructed.add(mt_str)
                added = True
                
            if not added:
                break
    
    return data

openmc.data.IncidentNeutron.add_hierarchical_mts = incident_neutron_add_hierarchical_mts

print("OpenMC monkey patches loaded successfully!")