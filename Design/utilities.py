from qiskit_metal import MetalGUI, Dict, Headings
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround

def location_str2float(pos):
    if pos[-2:] == 'um':
        return float(pos[:-2]) * 1e-3
    elif pos[-2:] == 'mm':
        return float(pos[:-2]) 
        
def convert_str_to_um(string):
    if string[-2:] == 'um':
        return float(string[:-2])

def add_length_in_str(alist):
    length = [convert_str_to_um(a) for a in alist]
    return str(sum(length)) + alist[0][-2:]


def connect(design, component_name: str, component1: str, pin1: str, component2: str, pin2: str,
            length: str, asymmetry='0 um', flip=False, fillet='99.99um',options=Dict()):
    """Connect two pins with a CPW."""
    myoptions = Dict(
        fillet=fillet,
        hfss_wire_bonds = True,
        pin_inputs=Dict(
            start_pin=Dict(
                component=component1,
                pin=pin1),
            end_pin=Dict(
                component=component2,
                pin=pin2)),
        total_length=length)

    myoptions.update(options)
    myoptions.meander.asymmetry = asymmetry
    myoptions.meander.lead_direction_inverted = 'true' if flip else 'false'
    return RouteMeander(design, component_name, myoptions)

