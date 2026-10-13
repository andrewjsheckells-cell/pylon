#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pylon
=====

A PySide6 tool for drawing a spacecraft liquid-propulsion feed-system schematic
and solving it as a steady, incompressible, adiabatic flow network.
(Called "Propulsion Schematic Builder" up to version 1.6.)

    pip install PySide6 numpy
    python pylon.py

See "Pylon.md" for the full user guide, the physics reference, the original
prompt and the living prompt.

Naming convention: identifiers are written out in full (no abbreviations).
Equations carry a comment naming their source from the reference list below.

References
----------
[CRANE]      Crane Co., "Flow of Fluids Through Valves, Fittings, and Pipe,"
             Technical Paper No. 410 (TP-410).  Darcy equation, resistance
             coefficient K, K = n * f_T fitting multipliers, pipe-bend table,
             sudden contraction / enlargement, entrance and exit losses,
             square-edged orifice discharge coefficient.
[CHURCHILL]  Churchill, S. W., "Friction-factor equation spans all
             fluid-flow regimes," Chemical Engineering 84(24), 91-92, 1977.
[COLEBROOK]  Colebrook, C. F., "Turbulent flow in pipes, with particular
             reference to the transition region between the smooth and rough
             pipe laws," J. Institution of Civil Engineers 11, 133-156, 1939.
[HOOPER]     Hooper, W. B., "The two-K method predicts head losses in pipe
             fittings," Chemical Engineering 88(17), 96-100, 1981.
[BLEVINS]    Blevins, R. D., Applied Fluid Dynamics Handbook, Van Nostrand
             Reinhold, 1984 (loss coefficients across laminar and turbulent flow).
[WHITE]      White, F. M., Fluid Mechanics, McGraw-Hill (Hagen-Poiseuille flow,
             Darcy-Weisbach equation, Borda-Carnot sudden-expansion loss,
             properties of water).
[SUTTON]     Sutton, G. P. and Biblarz, O., Rocket Propulsion Elements,
             9th ed., Wiley, 2017 (specific impulse and thrust; injector /
             orifice flow mdot = Cd * A * sqrt(2 * rho * dP)).
[ISA]        ANSI/ISA-75.01.01 (IEC 60534-2-1), Flow Equations for Sizing
             Control Valves (definition of the flow coefficient Cv).
[HYDRAZINE]  Arch Chemicals (Lonza), "Anhydrous Hydrazine / UPH / HPH Physical
             Properties," hydrazine.com Technical Data (density 1.004 g/mL and
             viscosity 0.90 cP at 25 C).
[NIST811]    Thompson, A. and Taylor, B. N., Guide for the Use of the
             International System of Units (SI), NIST Special Publication 811,
             2008 (unit conversion factors and standard gravity).
[TODINI]     Todini, E. and Pilati, S., "A gradient algorithm for the analysis
             of pipe networks," in Computer Applications in Water Supply,
             Vol. 1, Research Studies Press, 1988 (simultaneous solution of
             node heads and element flows).
[DENNIS]     Dennis, J. E. and Schnabel, R. B., Numerical Methods for
             Unconstrained Optimization and Nonlinear Equations, Prentice-Hall,
             1983 (Newton's method, backtracking line search, finite-difference
             derivatives).
[DOWELL]     Dowell, M. and Jarratt, P., "A modified regula falsi method for
             computing the root of an equation," BIT 11, 168-174, 1971
             (Illinois method).
[KIEFER]     Kiefer, J., "Sequential minimax search for a maximum,"
             Proceedings of the American Mathematical Society 4(3), 502-506,
             1953 (golden-section search).
[HART]       Hart, P. E., Nilsson, N. J. and Raphael, B., "A formal basis for
             the heuristic determination of minimum cost paths," IEEE Trans.
             Systems Science and Cybernetics 4(2), 100-107, 1968 (A* search).
[TARJAN]     Tarjan, R. E., "Efficiency of a good but not linear set union
             algorithm," Journal of the ACM 22(2), 215-225, 1975 (union-find
             with path compression).
"""
from __future__ import annotations

import copy
import csv
import heapq
import json
import math
import os
import re
import sys
import uuid

import numpy

APPLICATION_NAME = "Pylon"
APPLICATION_VERSION = "1.15"
FILE_FORMAT_VERSION = 3

# =============================================================================
# Units and constants  [NIST811]
# =============================================================================
STANDARD_GRAVITY = 9.80665                          # m/s^2, [NIST811]
PASCALS_PER_PSI = 6894.757293168                    # Pa per lbf/in^2, [NIST811]
METERS_PER_INCH = 0.0254                            # exact, [NIST811]
KILOGRAMS_PER_GRAM = 1.0e-3
CUBIC_METERS_PER_LITER = 1.0e-3
CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE = 1.0e-6 / 60.0
CUBIC_METERS_PER_SECOND_PER_GALLON_PER_MINUTE = 3.785411784e-3 / 60.0   # US gallon, [NIST811]
WATER_DENSITY_FOR_FLOW_COEFFICIENT = 999.0          # kg/m^3, water at 60 F used to define Cv, [ISA]

# Flow coefficient of a sharp-edged orifice, Cv = constant * Cd * d^2 (d in inches).
#   Orifice flow (volumetric):   Q = Cd * A * sqrt(2 * dP / rho)              [SUTTON], [CRANE]
#   Flow coefficient definition: Q[gpm] = Cv * sqrt(dP[psi] / specific gravity)   [ISA]
#   With dP = 1 psi of water:    Cv = Cd * (pi/4) d^2 * sqrt(2 * 1 psi / rho_water) / (1 gpm)
FLOW_COEFFICIENT_PER_DISCHARGE_COEFFICIENT_PER_SQUARE_INCH = (
    (math.pi / 4.0 * METERS_PER_INCH ** 2)
    * math.sqrt(2.0 * PASCALS_PER_PSI / WATER_DENSITY_FOR_FLOW_COEFFICIENT)
    / CUBIC_METERS_PER_SECOND_PER_GALLON_PER_MINUTE)       # = 29.84


def flow_coefficient_from_orifice_diameter(diameter_inches, discharge_coefficient):
    """Cv = 29.84 * Cd * d^2   [ISA] with orifice flow from [SUTTON]/[CRANE]."""
    return (FLOW_COEFFICIENT_PER_DISCHARGE_COEFFICIENT_PER_SQUARE_INCH
            * discharge_coefficient * diameter_inches * diameter_inches)


def orifice_diameter_from_flow_coefficient(flow_coefficient, discharge_coefficient):
    """Inverse of Cv = 29.84 * Cd * d^2   [ISA]."""
    if flow_coefficient <= 0 or discharge_coefficient <= 0:
        return 0.0
    return math.sqrt(flow_coefficient
                     / (FLOW_COEFFICIENT_PER_DISCHARGE_COEFFICIENT_PER_SQUARE_INCH * discharge_coefficient))


def format_number(value, significant_digits=4):
    """Compact number formatting for labels and tables."""
    if value is None:
        return "—"
    try:
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return "—"
    except TypeError:
        return str(value)
    if isinstance(value, int):
        return str(value)
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e5 or magnitude < 1e-3:
        return f"{value:.{significant_digits - 1}e}"
    return f"{value:.{significant_digits}g}"


# =============================================================================
# Tubing
# =============================================================================
TUBE_SIZES = [  # label, outer diameter [in], wall thickness [in]
    ('1/16" x 0.014"', 0.0625, 0.014),
    ('1/8" x 0.028"', 0.125, 0.028),
    ('1/8" x 0.020"', 0.125, 0.020),
    ('3/16" x 0.035"', 0.1875, 0.035),
    ('1/4" x 0.035"', 0.25, 0.035),
    ('1/4" x 0.049"', 0.25, 0.049),
    ('3/8" x 0.035"', 0.375, 0.035),
    ('3/8" x 0.049"', 0.375, 0.049),
    ('1/2" x 0.049"', 0.5, 0.049),
    ('1/2" x 0.065"', 0.5, 0.065),
    ('3/4" x 0.065"', 0.75, 0.065),
    ('1" x 0.083"', 1.0, 0.083),
]
TUBE_LABELS = [entry[0] for entry in TUBE_SIZES]
TUBE_DIMENSIONS_BY_LABEL = {entry[0]: (entry[1], entry[2]) for entry in TUBE_SIZES}
DEFAULT_TUBE = '1/4" x 0.035"'


def _parse_number(text):
    text = text.strip().strip('"').strip()
    if '/' in text:
        numerator, denominator = text.split('/', 1)
        return float(numerator) / float(denominator)
    return float(text)


def parse_tube(label):
    """(outer diameter, wall thickness) in inches for a standard label or a custom
    'outer diameter x wall' string, or None if the text cannot be understood."""
    if label in TUBE_DIMENSIONS_BY_LABEL:
        return TUBE_DIMENSIONS_BY_LABEL[label]
    match = re.match(r'\s*([\d./]+)\s*"?\s*[xX×]\s*([\d.]+)\s*"?\s*$', str(label))
    if match:
        try:
            outer_diameter = _parse_number(match.group(1))
            wall_thickness = _parse_number(match.group(2))
            if outer_diameter > 2 * wall_thickness > 0:
                return outer_diameter, wall_thickness
        except (ValueError, ZeroDivisionError):
            pass
    return None


def tube_outer_diameter_and_wall(label):
    return parse_tube(label) or TUBE_DIMENSIONS_BY_LABEL[DEFAULT_TUBE]


def tube_inner_diameter_inches(label):
    outer_diameter, wall_thickness = tube_outer_diameter_and_wall(label)
    return outer_diameter - 2.0 * wall_thickness


# =============================================================================
# Settings, fluids
# =============================================================================
FLUID_PRESETS = {
    # name: (density kg/m^3, dynamic viscosity mPa*s).
    # Hydrazine (the working fluid of this tool): anhydrous N2H4 at 25 C,
    #   1.004 g/mL and 0.90 cP                                              [HYDRAZINE]
    # Water from [WHITE] property tables; the other entries are approximate
    # handbook figures near 20 C, kept for ground-test fluids and comparisons.
    'Hydrazine (N2H4)': (1004.0, 0.90),
    'Monomethylhydrazine': (875.0, 0.86),
    'Nitrogen tetroxide (MON-3)': (1440.0, 0.42),
    'Water': (998.2, 1.002),
    'Isopropyl alcohol': (786.0, 2.40),
    'HFE-7100 (test fluid)': (1510.0, 0.58),
}
DEFAULT_FLUID = {'name': 'Hydrazine (N2H4)', 'density': 1004.0,             # [HYDRAZINE], 25 C
                 'viscosity_millipascal_seconds': 0.90}

SETTING_DEFINITIONS = [
    # key, label, default, tooltip
    ('orifice_discharge_coefficient', 'Orifice block discharge coefficient Cd', 0.6,
     'Cd used for Orifice blocks (square-edged orifice, [CRANE]; project instruction Cd = 0.6).'),
    ('crane_reference_roughness', 'Crane f_T roughness basis [in]', 0.0018,
     'Roughness used to evaluate the fully turbulent friction factor f_T that scales '
     'Crane K = n*f_T (Crane basis is clean commercial steel, 0.0018 in) [CRANE].'),
    ('tank_entrance_loss_coefficient', 'Tank entrance loss coefficient K (flow out of tank)', 0.5,
     'Square-edged entrance, K = 0.5 [CRANE].'),
    ('tank_exit_loss_coefficient', 'Tank exit loss coefficient K (flow into tank)', 1.0,
     'Exit into a large volume, K = 1.0 [CRANE].'),
    ('elbow_friction_multiplier', '90° fitting  K∞ = n·f_T,  n =', 30.0,
     'Crane standard 90° elbow: 30 f_T [CRANE].'),
    ('elbow_laminar_coefficient', '90° fitting laminar K1', 800.0,
     'Two-K laminar coefficient, K = K1/Re + K∞ [HOOPER].'),
    ('tee_run_friction_multiplier', 'Tee flow-through-run  n =', 20.0,
     'Crane tee, flow through run: 20 f_T [CRANE].'),
    ('tee_branch_friction_multiplier', 'Tee flow-through-branch  n =', 60.0,
     'Crane tee, flow through branch: 60 f_T [CRANE].'),
    ('tee_run_laminar_coefficient', 'Tee run laminar K1', 200.0,
     'Two-K laminar coefficient, run [HOOPER].'),
    ('tee_branch_laminar_coefficient', 'Tee branch laminar K1', 500.0,
     'Two-K laminar coefficient, branch [HOOPER].'),
    ('bend_laminar_coefficient', 'Tube bend laminar K1 (per bend)', 800.0,
     'Two-K laminar coefficient per 90° bend [HOOPER].'),
    ('size_mismatch_tolerance', 'Tube/port size mismatch tolerance [-]', 0.02,
     'Relative inner-diameter difference above which contraction/expansion losses are added.'),
    ('acceptance_tolerance_percent', 'Thruster flow acceptance tolerance [%]', 2.0,
     'Used by orifice sizing and the lowest-RPM search.'),
]
DEFAULT_SETTINGS = {key: default for key, _, default, _ in SETTING_DEFINITIONS}

# 90° pipe bends: K = n * f_T as a function of bend radius / diameter   [CRANE]
CRANE_BEND_MULTIPLIER_BY_RADIUS_RATIO = [(1, 20), (1.5, 14), (2, 12), (3, 12), (4, 14), (6, 17),
                                         (8, 24), (10, 30), (12, 34), (14, 38), (16, 42), (20, 50)]


def crane_bend_multiplier(bend_radius_ratio):
    """n in K = n * f_T for a 90° bend, linearly interpolated in the [CRANE] table."""
    radius_ratios = [ratio for ratio, _ in CRANE_BEND_MULTIPLIER_BY_RADIUS_RATIO]
    multipliers = [multiplier for _, multiplier in CRANE_BEND_MULTIPLIER_BY_RADIUS_RATIO]
    return float(numpy.interp(bend_radius_ratio, radius_ratios, multipliers))


def churchill_friction_factor_times_reynolds(reynolds_number, relative_roughness):
    """Darcy friction factor multiplied by the Reynolds number   [CHURCHILL].

        f = 8 * [ (8/Re)^12 + (A + B)^(-3/2) ]^(1/12)
        A = [ -2.457 * ln( (7/Re)^0.9 + 0.27 * epsilon/D ) ]^16
        B = (37530/Re)^16

    f*Re tends to 64 in laminar flow (Hagen-Poiseuille, [WHITE]), so the
    friction pressure drop written as (f*Re) * mu * V * L / (2 D^2) is smooth
    through zero flow.
    """
    if reynolds_number < 1.0:
        return 64.0                                   # laminar limit f = 64/Re   [WHITE]
    laminar_term = (8.0 / reynolds_number) ** 12
    term_a = (-2.457 * math.log((7.0 / reynolds_number) ** 0.9 + 0.27 * relative_roughness)) ** 16
    term_b = (37530.0 / reynolds_number) ** 16
    friction_factor = 8.0 * (laminar_term + (term_a + term_b) ** -1.5) ** (1.0 / 12.0)   # [CHURCHILL]
    return friction_factor * reynolds_number


def churchill_friction_factor(reynolds_number, relative_roughness):
    """Darcy friction factor f   [CHURCHILL]."""
    return (churchill_friction_factor_times_reynolds(reynolds_number, relative_roughness)
            / max(reynolds_number, 1e-300))


def fully_turbulent_friction_factor(diameter_meters, settings):
    """f_T used to scale Crane fitting coefficients K = n * f_T   [CRANE].

    Evaluated as the Re -> infinity limit of [CHURCHILL] (B -> 0, (8/Re)^12 -> 0):
        f_T = 8 * A^(-1/8) = 8 / [2.457 * ln(0.27 * epsilon/D)]^2
    on the Crane roughness basis (clean commercial steel)."""
    relative_roughness = (settings.get('crane_reference_roughness', 0.0018) * METERS_PER_INCH
                          / max(diameter_meters, 1e-9))
    if relative_roughness <= 0:
        return 0.008
    return 1.0 / (2.457 / (2.0 * math.sqrt(2.0)) * math.log(0.27 * relative_roughness)) ** 2


# =============================================================================
# Block / line property definitions
# =============================================================================
class PropertyDefinition:
    """Definition of one editable or calculated property of a block or line."""
    __slots__ = ('key', 'label', 'unit', 'default', 'kind', 'calculated', 'choices', 'tooltip')

    def __init__(self, key, label, unit='', default=0.0, kind='float', calculated=False,
                 choices=None, tooltip=''):
        self.key, self.label, self.unit, self.default = key, label, unit, default
        self.kind, self.calculated, self.choices, self.tooltip = kind, calculated, choices, tooltip


def _tube_property(key, label, default=DEFAULT_TUBE):
    return PropertyDefinition(key, label, '', default, 'tube',
                              tooltip='Standard size or custom "outer diameter x wall" in inches.')


def _calculated_property(key, label, unit=''):
    return PropertyDefinition(key, label, unit, None, 'float', calculated=True)


BLOCK_TYPES = {
    'tank': dict(
        # drawn 4x the earlier tank: 176 x 280 capsule plus inlet / outlet tube stubs
        label='Tank', prefix='TK', size=(200, 400), color='#3d6fb6',
        ports={'inlet': (0, -200, 0, -1), 'outlet': (0, 200, 0, 1)},
        states=['Nominal pressure', 'Maximum pressure', 'Minimum pressure'],
        properties=[PropertyDefinition('volume', 'Volume', 'L', 10.0),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure', 'psia', 60.0),
                    PropertyDefinition('nominal_pressure', 'Nominal pressure', 'psia', 50.0),
                    PropertyDefinition('minimum_pressure', 'Minimum pressure', 'psia', 30.0),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size (top)'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size (bottom)'),
                    _calculated_property('calculated_pressure', 'Tank pressure (current state)', 'psia'),
                    _calculated_property('calculated_inflow', 'Inflow (top inlet)', 'g/s'),
                    _calculated_property('calculated_net_outflow', 'Net outflow', 'g/s'),
                    _calculated_property('calculated_time_to_empty', 'Time to empty at current outflow', 's')],
        flow_coefficient_links=[]),
    'copv': dict(
        # composite overwrapped pressure vessel: pressurant bottle that backfills tank ullage
        label='COPV', prefix='COPV', size=(100, 140),
        symbol_extent=(55, 70, 30, 30), color='#8e7cc3',   # extent: symbol reach up/down/left/right
        ports={'outlet': (0, 70, 0, 1)},
        states=['Backfilling', 'Isolated'],
        properties=[PropertyDefinition('volume', 'Volume', 'L', 10.0),
                    PropertyDefinition('bottle_pressure', 'Bottle pressure', 'psia', 3000.0),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure (MEOP)', 'psia', 4500.0),
                    PropertyDefinition('ullage_pressure_setpoint', 'Ullage pressure setpoint (regulated)', 'psia',
                                       50.0, tooltip='Pressure the COPV holds in the tank ullage. '
                                                     '0 = unregulated (full bottle pressure).'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    _calculated_property('calculated_delivered_pressure', 'Ullage pressure delivered', 'psia'),
                    _calculated_property('calculated_tanks_backfilled', 'Tanks backfilled', ''),
                    _calculated_property('calculated_backfill_volume_flow', 'Ullage backfill volume flow',
                                         'mL/min')],
        flow_coefficient_links=[]),
    'latch_valve': dict(
        label='Latch Valve', prefix='LV', size=(120, 100),
        symbol_extent=(25, 13, 22, 22), color='#7a4fb0',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['Open', 'Closed'],
        properties=[PropertyDefinition('open_orifice_diameter', 'Equivalent orifice diameter (open state)', 'in', 0.10),
                    PropertyDefinition('closed_orifice_diameter', 'Equivalent orifice diameter (closed state, leak)',
                                       'in', 0.0, tooltip='0 = bubble-tight'),
                    PropertyDefinition('discharge_coefficient', 'Discharge coefficient', '', 0.6),
                    PropertyDefinition('flow_coefficient', 'Flow coefficient Cv (open)', 'gpm/√psi', 0.0,
                                       tooltip='Linked to the equivalent orifice diameter and discharge coefficient'),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing diameter'),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure', 'psia', 400.0),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_drop', 'Pressure drop', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Mass flow', 'g/s')],
        flow_coefficient_links=[('open_orifice_diameter', 'flow_coefficient', 'discharge_coefficient')]),
    'fluid_coupling': dict(
        # quick-disconnect style coupling: open (connected) with a flow restriction, or closed
        label='Fluid Coupling', prefix='FC', size=(120, 90),
        symbol_extent=(9, 9, 24, 24), color='#3a7d9c',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['Open', 'Closed'],
        properties=[PropertyDefinition('orifice_diameter', 'Equivalent orifice diameter (open, flow restriction)',
                                       'in', 0.125),
                    PropertyDefinition('discharge_coefficient', 'Discharge coefficient', '', 0.6),
                    PropertyDefinition('flow_coefficient', 'Flow coefficient Cv (open)', 'gpm/√psi', 0.0,
                                       tooltip='Linked to the equivalent orifice diameter and discharge coefficient'),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure', 'psia', 400.0),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_drop', 'Pressure drop', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Mass flow', 'g/s')],
        flow_coefficient_links=[('orifice_diameter', 'flow_coefficient', 'discharge_coefficient')]),
    'fill_drain_valve': dict(
        # symbolic service valve: connects like a component but passes no flow (dead end each side)
        label='Fill & Drain Valve', prefix='FDV', size=(120, 90),
        symbol_extent=(18, 13, 22, 22), color='#7f8c8d',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['Closed (dead end)'],
        properties=[_tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure', 'psia', 400.0),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia')],
        flow_coefficient_links=[]),
    'filter': dict(
        label='Filter', prefix='FL', size=(120, 110),
        symbol_extent=(19, 19, 20, 20), color='#2f8f6a',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['Clean', 'Maximum trapped particles'],
        properties=[PropertyDefinition('clean_orifice_diameter', 'Equivalent orifice diameter (no trapped particles)',
                                       'in', 0.06),
                    PropertyDefinition('clogged_orifice_diameter',
                                       'Equivalent orifice diameter (maximum trapped particles)', 'in', 0.035),
                    PropertyDefinition('discharge_coefficient', 'Discharge coefficient', '', 0.6),
                    PropertyDefinition('clean_flow_coefficient', 'Flow coefficient Cv (clean)', 'gpm/√psi', 0.0),
                    PropertyDefinition('clogged_flow_coefficient', 'Flow coefficient Cv (maximum trapped)',
                                       'gpm/√psi', 0.0),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure', 'psia', 400.0),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_drop', 'Pressure drop', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Mass flow', 'g/s')],
        flow_coefficient_links=[('clean_orifice_diameter', 'clean_flow_coefficient', 'discharge_coefficient'),
                                ('clogged_orifice_diameter', 'clogged_flow_coefficient', 'discharge_coefficient')]),
    'flowmeter': dict(
        label='Flowmeter', prefix='FM', size=(120, 100),
        symbol_extent=(17, 17, 17, 17), color='#b0742a',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['Flowing', 'Blocked'],
        properties=[PropertyDefinition('orifice_diameter', 'Equivalent orifice diameter', 'in', 0.07),
                    PropertyDefinition('flow_coefficient', 'Flow coefficient Cv', 'gpm/√psi', 0.0),
                    PropertyDefinition('discharge_coefficient', 'Discharge coefficient', '', 0.6),
                    PropertyDefinition('maximum_pressure', 'Maximum pressure', 'psia', 400.0),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_drop', 'Pressure drop', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Measured mass flow', 'g/s')],
        flow_coefficient_links=[('orifice_diameter', 'flow_coefficient', 'discharge_coefficient')]),
    'transducer': dict(
        label='Pressure Transducer', prefix='PT', size=(80, 80),
        symbol_extent=(25, 9, 17, 17), color='#c0392b',   # extent: symbol reach up/down/left/right
        ports={'inlet': (0, 40, 0, 1)},
        states=['ON', 'OFF'],
        properties=[_tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _calculated_property('calculated_pressure', 'Pressure reading', 'psia')],
        flow_coefficient_links=[]),
    'pump': dict(
        label='Pump', prefix='P', size=(120, 120),
        symbol_extent=(21, 21, 21, 21), color='#16758f',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['ON', 'OFF'],
        properties=[PropertyDefinition('revolutions_per_minute', 'Speed', 'rpm', 3000.0,
                                       tooltip='When ON the pump follows its imported curve for this speed '
                                               '(interpolated between stored speeds); also editable on the block.'),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    PropertyDefinition('stopped_orifice_diameter', 'Equivalent orifice diameter (stopped)', 'in', 0.02),
                    PropertyDefinition('stopped_discharge_coefficient', 'Discharge coefficient (stopped)', '', 0.6),
                    PropertyDefinition('stopped_flow_coefficient', 'Flow coefficient Cv (stopped)', 'gpm/√psi', 0.0),
                    PropertyDefinition('maximum_inlet_pressure', 'Maximum inlet pressure', 'psia', 100.0),
                    PropertyDefinition('minimum_inlet_pressure', 'Minimum inlet pressure', 'psia', 15.0),
                    PropertyDefinition('maximum_outlet_pressure', 'Maximum outlet pressure', 'psia', 450.0),
                    PropertyDefinition('maximum_differential_pressure', 'Maximum differential pressure', 'psid', 400.0),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_rise', 'Pressure rise', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Mass flow', 'g/s'),
                    _calculated_property('calculated_volumetric_flow', 'Volumetric flow', 'mL/min')],
        flow_coefficient_links=[('stopped_orifice_diameter', 'stopped_flow_coefficient',
                                 'stopped_discharge_coefficient')]),
    'thruster': dict(
        label='Thruster', prefix='TH', size=(120, 110),
        symbol_extent=(21, 21, 34, 52), color='#d35400',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0)},
        states=['ON', 'OFF'],
        properties=[PropertyDefinition('specific_impulse', 'Specific impulse', 's', 230.0),
                    PropertyDefinition('nominal_inlet_pressure', 'Nominal inlet pressure', 'psia', 300.0),
                    PropertyDefinition('maximum_inlet_pressure', 'Maximum inlet pressure', 'psia', 400.0),
                    PropertyDefinition('orifice_diameter', 'Equivalent orifice diameter', 'in', 0.0143),
                    PropertyDefinition('discharge_coefficient', 'Discharge coefficient', '', 0.6),
                    PropertyDefinition('flow_coefficient', 'Flow coefficient Cv', 'gpm/√psi', 0.0),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size', '1/8" x 0.028"'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure (hard coded)', 'psia'),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_mass_flow', 'Flowrate', 'g/s'),
                    _calculated_property('calculated_thrust', 'Thrust', 'N'),
                    _calculated_property('calculated_nominal_mass_flow', 'Flowrate at nominal inlet pressure', 'g/s'),
                    _calculated_property('calculated_nominal_thrust', 'Thrust at nominal inlet pressure', 'N')],
        flow_coefficient_links=[('orifice_diameter', 'flow_coefficient', 'discharge_coefficient')]),
    'orifice': dict(
        label='Orifice', prefix='OR', size=(120, 110),
        symbol_extent=(20, 20, 4, 4), color='#566573',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-60, 0, -1, 0), 'outlet': (60, 0, 1, 0)},
        states=['ON', 'OFF'],
        properties=[PropertyDefinition('orifice_diameter', 'Orifice diameter', 'in', 0.010),
                    _tube_property('inlet_tubing_size', 'Inlet tubing size'),
                    _tube_property('outlet_tubing_size', 'Outlet tubing size'),
                    _calculated_property('calculated_discharge_coefficient', 'Discharge coefficient (fixed)', ''),
                    _calculated_property('calculated_flow_coefficient', 'Flow coefficient Cv', 'gpm/√psi'),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_drop', 'Pressure drop', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Mass flow (inlet→outlet)', 'g/s')],
        flow_coefficient_links=[]),
    'tee': dict(
        label='Tee', prefix='T', size=(80, 80),
        symbol_extent=(4, 40, 40, 40), color='#5d6d7e',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-40, 0, -1, 0), 'outlet': (40, 0, 1, 0), 'branch': (0, 40, 0, 1)},
        states=['Open'],                       # fittings are always open: no toggle
        properties=[_tube_property('inlet_tubing_size', 'Inlet tube size'),
                    _tube_property('outlet_tubing_size', 'Outlet tube size'),
                    _tube_property('branch_tubing_size', 'Tee (branch) tube size'),
                    _calculated_property('calculated_center_pressure', 'Centre pressure', 'psia'),
                    _calculated_property('calculated_inlet_flow', 'Flow in (inlet port)', 'g/s'),
                    _calculated_property('calculated_outlet_flow', 'Flow out (outlet port)', 'g/s'),
                    _calculated_property('calculated_branch_flow', 'Flow out (branch port)', 'g/s')],
        flow_coefficient_links=[]),
    'cross': dict(
        # four-way fitting: run inlet -> outlet, branches up and down
        label='Cross', prefix='X', size=(80, 80),
        symbol_extent=(40, 40, 40, 40), color='#5d6d7e',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-40, 0, -1, 0), 'outlet': (40, 0, 1, 0),
               'upper_branch': (0, -40, 0, -1), 'lower_branch': (0, 40, 0, 1)},
        states=['Open'],                       # fittings are always open: no toggle
        properties=[_tube_property('inlet_tubing_size', 'Inlet tube size'),
                    _tube_property('outlet_tubing_size', 'Outlet tube size'),
                    _tube_property('upper_branch_tubing_size', 'Upper branch tube size'),
                    _tube_property('lower_branch_tubing_size', 'Lower branch tube size'),
                    _calculated_property('calculated_center_pressure', 'Centre pressure', 'psia'),
                    _calculated_property('calculated_inlet_flow', 'Flow in (inlet port)', 'g/s'),
                    _calculated_property('calculated_outlet_flow', 'Flow out (outlet port)', 'g/s'),
                    _calculated_property('calculated_upper_branch_flow', 'Flow out (upper branch)', 'g/s'),
                    _calculated_property('calculated_lower_branch_flow', 'Flow out (lower branch)', 'g/s')],
        flow_coefficient_links=[]),
    'elbow': dict(
        label='90° Fitting', prefix='EL', size=(80, 80),
        symbol_extent=(4, 40, 40, 4), color='#5d6d7e',   # extent: symbol reach up/down/left/right
        ports={'inlet': (-40, 0, -1, 0), 'outlet': (0, 40, 0, 1)},
        states=['Open'],                       # fittings are always open: no toggle
        properties=[_tube_property('inlet_tubing_size', 'Inlet tube size'),
                    _tube_property('outlet_tubing_size', 'Outlet tube size'),
                    _calculated_property('calculated_inlet_pressure', 'Inlet pressure', 'psia'),
                    _calculated_property('calculated_outlet_pressure', 'Outlet pressure', 'psia'),
                    _calculated_property('calculated_pressure_drop', 'Pressure drop', 'psid'),
                    _calculated_property('calculated_mass_flow', 'Mass flow', 'g/s')],
        flow_coefficient_links=[]),
}
BLOCK_ORDER = ['tank', 'copv', 'latch_valve', 'fluid_coupling', 'filter', 'flowmeter', 'fill_drain_valve',
               'transducer', 'pump',
               'thruster', 'orifice', 'tee', 'cross', 'elbow']

PORT_TUBING_PROPERTY = {
    'tank': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'copv': {'outlet': 'outlet_tubing_size'},
    'fluid_coupling': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'fill_drain_valve': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'latch_valve': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'filter': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'flowmeter': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'transducer': {'inlet': 'inlet_tubing_size'},
    'pump': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'thruster': {'inlet': 'inlet_tubing_size'},
    'orifice': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
    'tee': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size', 'branch': 'branch_tubing_size'},
    'cross': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size',
              'upper_branch': 'upper_branch_tubing_size', 'lower_branch': 'lower_branch_tubing_size'},
    'elbow': {'inlet': 'inlet_tubing_size', 'outlet': 'outlet_tubing_size'},
}

LINE_PROPERTIES = [
    PropertyDefinition('outer_diameter', 'Outer diameter', 'in', 0.25),
    PropertyDefinition('inner_diameter', 'Inner diameter', 'in', 0.18),
    PropertyDefinition('length', 'Length', 'in', 12.0),
    PropertyDefinition('bend_count', 'Number of bends (90°)', '', 0, 'int',
                       tooltip='Used in the flow calculation only; not tied to the drawn shape.'),
    PropertyDefinition('bend_radius_ratio', 'Bend radius / diameter', '', 3.0),
    PropertyDefinition('roughness', 'Absolute roughness', 'in', 0.00006,
                       tooltip='Drawn tubing is about 6e-5 in (1.5 µm).'),
    _calculated_property('calculated_mass_flow', 'Mass flow (start→end)', 'g/s'),
    _calculated_property('calculated_pressure_drop', 'Pressure drop inlet→outlet', 'psid'),
    _calculated_property('calculated_velocity', 'Velocity', 'm/s'),
    _calculated_property('calculated_reynolds_number', 'Reynolds number', ''),
]


def default_block_properties(block_type):
    properties = {definition.key: definition.default
                  for definition in BLOCK_TYPES[block_type]['properties'] if not definition.calculated}
    update_linked_flow_coefficients(block_type, properties)
    return properties


def default_line_properties():
    return {definition.key: definition.default for definition in LINE_PROPERTIES if not definition.calculated}


def update_linked_flow_coefficients(block_type, properties, changed_key=None):
    """Keep orifice diameter, discharge coefficient and Cv consistent
    (Cv = 29.84 * Cd * d^2, [ISA]).  If `changed_key` is a Cv key, the diameter follows."""
    for diameter_key, flow_coefficient_key, discharge_key in BLOCK_TYPES[block_type]['flow_coefficient_links']:
        discharge_coefficient = float(properties.get(discharge_key, 0.6) or 0.6)
        if changed_key == flow_coefficient_key:
            properties[diameter_key] = orifice_diameter_from_flow_coefficient(
                float(properties.get(flow_coefficient_key, 0.0)), discharge_coefficient)
        else:
            properties[flow_coefficient_key] = flow_coefficient_from_orifice_diameter(
                float(properties.get(diameter_key, 0.0)), discharge_coefficient)
    return properties


def state_names(block_type):
    return BLOCK_TYPES[block_type]['states']


def has_toggle(block_type):
    """True when the block has more than one state (tees and 90° fittings have one)."""
    return len(BLOCK_TYPES[block_type]['states']) > 1


def normalized_state(block_type, state):
    """Clamp a stored state index to the states this block type has."""
    state_count = len(BLOCK_TYPES[block_type]['states'])
    state = int(state or 0)
    return state if 0 <= state < state_count else 0


def is_on(block):
    """The first state of every block type is its normal / active one."""
    return int(block.get('state', 0)) == 0


TANK_PRESSURE_KEY_BY_STATE = ['nominal_pressure', 'maximum_pressure', 'minimum_pressure']


def tank_pressure_psia(properties, state):
    return float(properties.get(TANK_PRESSURE_KEY_BY_STATE[int(state) % 3], 0.0))


def copv_delivered_pressure(properties, state):
    """Ullage pressure a COPV delivers [psia], or None while isolated.  An ideal pressure
    regulator holds the setpoint (regulated pressure feed, [SUTTON]); it cannot deliver more
    than the bottle pressure.  A setpoint of 0 means unregulated (bottle pressure)."""
    if int(state or 0) != 0:
        return None
    bottle_pressure = float(properties.get('bottle_pressure', 0.0))
    setpoint = float(properties.get('ullage_pressure_setpoint', 0.0) or 0.0)
    return min(setpoint, bottle_pressure) if setpoint > 0 else bottle_pressure


def passes_pressurant(block):
    """True if pressurant can pass through this block between a COPV and a tank ullage:
    fittings and filters always, valves / couplings / flowmeters / orifices only when open."""
    block_type = block.get('type')
    state = normalized_state(block_type, block.get('state', 0)) if block_type in BLOCK_TYPES else 0
    if block_type in ('tee', 'cross', 'elbow', 'filter'):
        return True
    return block_type in ('latch_valve', 'fluid_coupling', 'flowmeter', 'orifice') and state == 0


def ullage_backfill(model):
    """{tank identifier: (ullage pressure psia, COPV identifier)} for every tank whose top
    inlet is connected to a backfilling COPV through lines and open passive components.
    With several COPVs on one tank the highest delivered pressure wins."""
    blocks = blocks_by_identifier(model)
    copvs = [block for block in blocks.values() if block.get('type') == 'copv']
    if not copvs:
        return {}
    connections = UnionFind()
    for line in model.get('lines', []):
        line_node = f"line:{line['identifier']}"
        connections.find(line_node)
        for end_name in ('start', 'end'):
            attachment = line.get(end_name) or {}
            if attachment.get('kind') == 'port' and attachment.get('block') in blocks:
                connections.union(line_node, f"{attachment['block']}:{attachment['port']}")
            elif attachment.get('kind') == 'line':
                connections.union(line_node, f"line:{attachment.get('line')}")
    for block in blocks.values():
        if block.get('type') in BLOCK_TYPES and passes_pressurant(block):
            port_nodes = [f"{block['identifier']}:{port}" for port in BLOCK_TYPES[block['type']]['ports']]
            for port_node in port_nodes[1:]:
                connections.union(port_nodes[0], port_node)
    backfill = {}
    for copv in copvs:
        delivered = copv_delivered_pressure(copv.get('properties', {}), copv.get('state', 0))
        if delivered is None:
            continue
        copv_root = connections.find(f"{copv['identifier']}:outlet")
        for tank in blocks.values():
            if tank.get('type') == 'tank' and connections.find(f"{tank['identifier']}:inlet") == copv_root:
                if tank['identifier'] not in backfill or delivered > backfill[tank['identifier']][0]:
                    backfill[tank['identifier']] = (delivered, copv['identifier'])
    return backfill


# =============================================================================
# Older file formats (versions 1 and 2 used abbreviated keys)
# =============================================================================
COMMON_PROPERTY_RENAMES = {'in_tube': 'inlet_tubing_size', 'out_tube': 'outlet_tubing_size',
                           'tee_tube': 'branch_tubing_size', 'p_max': 'maximum_pressure',
                           'cd': 'discharge_coefficient', 'cv': 'flow_coefficient'}
PROPERTY_RENAMES_BY_TYPE = {
    'tank': {'p_nom': 'nominal_pressure', 'p_min': 'minimum_pressure'},
    'latch_valve': {'eseod_open': 'open_orifice_diameter', 'eseod_closed': 'closed_orifice_diameter'},
    'filter': {'eseod': 'clean_orifice_diameter', 'eseod_dirty': 'clogged_orifice_diameter',
               'cv': 'clean_flow_coefficient', 'cv_dirty': 'clogged_flow_coefficient'},
    'flowmeter': {'eseod': 'orifice_diameter'},
    'transducer': {},
    'pump': {'mode': 'operating_model', 'flow': 'flow_rate', 'rpm': 'revolutions_per_minute',
             'eseod_stop': 'stopped_orifice_diameter', 'cd_stop': 'stopped_discharge_coefficient',
             'cv_stop': 'stopped_flow_coefficient', 'p_in_max': 'maximum_inlet_pressure',
             'p_in_min': 'minimum_inlet_pressure', 'p_out_max': 'maximum_outlet_pressure',
             'dp_max': 'maximum_differential_pressure'},
    'thruster': {'isp': 'specific_impulse', 'p_nom': 'nominal_inlet_pressure',
                 'p_max': 'maximum_inlet_pressure', 'eseod': 'orifice_diameter'},
    'orifice': {'diameter': 'orifice_diameter'},
    'tee': {},
    'cross': {},
    'elbow': {},
}
LINE_PROPERTY_RENAMES = {'od': 'outer_diameter', 'id': 'inner_diameter', 'bends': 'bend_count',
                         'bend_rd': 'bend_radius_ratio', 'rough': 'roughness'}
SETTING_RENAMES = {'orifice_cd': 'orifice_discharge_coefficient', 'ft_rough_in': 'crane_reference_roughness',
                   'tank_entr_K': 'tank_entrance_loss_coefficient', 'tank_exit_K': 'tank_exit_loss_coefficient',
                   'elbow_n': 'elbow_friction_multiplier', 'elbow_K1': 'elbow_laminar_coefficient',
                   'tee_run_n': 'tee_run_friction_multiplier', 'tee_branch_n': 'tee_branch_friction_multiplier',
                   'tee_run_K1': 'tee_run_laminar_coefficient', 'tee_branch_K1': 'tee_branch_laminar_coefficient',
                   'bend_K1': 'bend_laminar_coefficient', 'mismatch_tol': 'size_mismatch_tolerance',
                   'accept_tol_pct': 'acceptance_tolerance_percent'}
PORT_RENAMES = {'in': 'inlet', 'out': 'outlet'}


# Properties that no longer exist (removed from older versions); dropped when a file is read.
#   pump: fixed-flow mode was removed in v1.9 - an ON pump always follows its RPM curves.
OBSOLETE_PROPERTIES_BY_TYPE = {'pump': ('operating_model', 'flow_rate')}


def remove_obsolete_properties(block_type, properties):
    for key in OBSOLETE_PROPERTIES_BY_TYPE.get(block_type, ()):
        (properties or {}).pop(key, None)
    return properties


def remove_obsolete_model_properties(model):
    """Drop obsolete properties from a current-version model, in place."""
    type_by_identifier = {}
    for block in model.get('blocks', []):
        remove_obsolete_properties(block.get('type'), block.get('properties'))
        type_by_identifier[block.get('identifier')] = block.get('type')
    for scenario in model.get('scenarios', []) or []:
        for identifier, snapshot in (scenario.get('blocks') or {}).items():
            remove_obsolete_properties(type_by_identifier.get(identifier), snapshot.get('properties'))
    return model


def _rename_keys(dictionary, renames):
    return {renames.get(key, key): value for key, value in (dictionary or {}).items()}


def upgrade_block_properties(block_type, properties):
    renames = dict(COMMON_PROPERTY_RENAMES)
    renames.update(PROPERTY_RENAMES_BY_TYPE.get(block_type, {}))
    upgraded = _rename_keys(properties, renames)
    upgraded = {key: value for key, value in upgraded.items() if not key.startswith('c_')}
    return remove_obsolete_properties(block_type, upgraded)


def upgrade_block_template(template):
    """Saved-block JSON from any version -> current keys."""
    upgraded = dict(template)
    if 'props' in upgraded:
        upgraded['properties'] = upgrade_block_properties(upgraded['type'], upgraded.pop('props'))
    if 'curves' in upgraded:
        upgraded['pump_curves'] = upgraded.pop('curves')
    if 'orient' in upgraded:
        upgraded['orientation'] = upgraded.pop('orient')
    upgraded['properties'] = remove_obsolete_properties(upgraded.get('type'),
                                                        dict(upgraded.get('properties') or {}))
    upgraded['format'] = 'pylon_block'
    upgraded['version'] = FILE_FORMAT_VERSION
    return upgraded


def _upgrade_attachment(attachment):
    if not attachment:
        return attachment
    upgraded = dict(attachment)
    if 'port' in upgraded:
        upgraded['port'] = PORT_RENAMES.get(upgraded['port'], upgraded['port'])
    if 't' in upgraded:
        upgraded['path_fraction'] = upgraded.pop('t')
    return upgraded


def upgrade_model(model):
    """Schematic JSON from any earlier version -> the current format (version 3)."""
    if int(model.get('version', 1)) >= FILE_FORMAT_VERSION:
        return remove_obsolete_model_properties(model)
    upgraded = copy.deepcopy(model)
    fluid = upgraded.get('fluid') or {}
    upgraded['fluid'] = {'name': fluid.get('name', DEFAULT_FLUID['name']),
                         'density': fluid.get('rho', fluid.get('density', DEFAULT_FLUID['density'])),
                         'viscosity_millipascal_seconds': fluid.get(
                             'mu_mPas', fluid.get('viscosity_millipascal_seconds',
                                                  DEFAULT_FLUID['viscosity_millipascal_seconds']))}
    upgraded['settings'] = _rename_keys(upgraded.get('settings'), SETTING_RENAMES)
    type_by_identifier = {}
    for block in upgraded.get('blocks', []):
        if 'id' in block:
            block['identifier'] = block.pop('id')
        if 'orient' in block:
            block['orientation'] = block.pop('orient')
        if 'curves' in block:
            block['pump_curves'] = block.pop('curves')
        block['properties'] = upgrade_block_properties(block['type'], block.pop('props', block.get('properties')))
        type_by_identifier[block['identifier']] = block['type']
    for line in upgraded.get('lines', []):
        if 'id' in line:
            line['identifier'] = line.pop('id')
        if 'auto' in line:
            line['auto_routed'] = line.pop('auto')
        line['properties'] = _rename_keys(line.pop('props', line.get('properties')), LINE_PROPERTY_RENAMES)
        line['start'] = _upgrade_attachment(line.get('start'))
        line['end'] = _upgrade_attachment(line.get('end'))
    for scenario in upgraded.get('scenarios', []):
        for identifier, snapshot in scenario.get('blocks', {}).items():
            snapshot['properties'] = upgrade_block_properties(type_by_identifier.get(identifier, ''),
                                                              snapshot.pop('props', snapshot.get('properties')))
        for snapshot in scenario.get('lines', {}).values():
            snapshot['properties'] = _rename_keys(snapshot.pop('props', snapshot.get('properties')),
                                                  LINE_PROPERTY_RENAMES)
    upgraded['format'] = 'pylon_schematic'
    upgraded['version'] = FILE_FORMAT_VERSION
    return upgraded


# =============================================================================
# Pump curves
# =============================================================================
def curve_points_for_speed(pump_curves, revolutions_per_minute):
    """Curve stored for exactly this speed (to the nearest rpm), else None."""
    return pump_curves.get(str(int(round(revolutions_per_minute)))) if pump_curves else None


def make_pump_curve_function_for_speed(pump_curves, revolutions_per_minute, density):
    """Pump curve at any speed.  Returns (function, note) where function is as for
    make_pump_curve_function and note is a warning text or None.

    * A stored speed uses its own curve.
    * A speed between two stored speeds N_low < N < N_high interpolates linearly in
      speed at the same pressure rise:
          Q(dP, N) = (1 - w) * Q_low(dP) + w * Q_high(dP),   w = (N - N_low) / (N_high - N_low)
      (exact for a positive-displacement pump, Q = displacement * N - slip(dP)).
    * A speed outside the stored range uses the nearest end curve and returns a warning.
    Before v1.6 an in-between speed silently used the nearest stored curve."""
    if not pump_curves:
        return None, None
    exact_points = curve_points_for_speed(pump_curves, revolutions_per_minute)
    if exact_points:
        return make_pump_curve_function(exact_points, density), None
    speeds = sorted((float(key), key) for key in pump_curves if pump_curves[key])
    if not speeds:
        return None, None
    if revolutions_per_minute <= speeds[0][0] or revolutions_per_minute >= speeds[-1][0]:
        end_speed, end_key = speeds[0] if revolutions_per_minute <= speeds[0][0] else speeds[-1]
        note = (f"Speed {revolutions_per_minute:g} rpm is outside the imported curves "
                f"({speeds[0][0]:g}–{speeds[-1][0]:g} rpm); using the {end_speed:g} rpm curve.")
        return make_pump_curve_function(pump_curves[end_key], density), note
    for (low_speed, low_key), (high_speed, high_key) in zip(speeds, speeds[1:]):
        if low_speed < revolutions_per_minute < high_speed:
            weight = (revolutions_per_minute - low_speed) / (high_speed - low_speed)
            low_curve = make_pump_curve_function(pump_curves[low_key], density)
            high_curve = make_pump_curve_function(pump_curves[high_key], density)

            def interpolated_curve(pressure_rise, low_curve=low_curve, high_curve=high_curve, weight=weight):
                low_flow, low_slope = low_curve(pressure_rise)
                high_flow, high_slope = high_curve(pressure_rise)
                return ((1 - weight) * low_flow + weight * high_flow,
                        (1 - weight) * low_slope + weight * high_slope)
            return interpolated_curve, None
    return None, None


def make_pump_curve_function(points, density):
    """Return f(pressure_rise_pascals) -> (mass_flow kg/s, d mass_flow / d pressure_rise)
    from [(pressure_rise_psid, volumetric_flow_mL_per_min), ...].

    Piecewise-linear interpolation, linearly extrapolated past the end points;
    mass flow = density * volumetric flow (incompressible fluid)."""
    sorted_points = sorted((float(pressure_rise), float(flow)) for pressure_rise, flow in points)
    if len(sorted_points) == 1:
        sorted_points = [sorted_points[0], (sorted_points[0][0] + 1.0, sorted_points[0][1])]
    pressure_rises = numpy.array([point[0] * PASCALS_PER_PSI for point in sorted_points])
    mass_flows = numpy.array([point[1] * CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE * density
                              for point in sorted_points])

    def pump_curve(pressure_rise):
        if pressure_rise <= pressure_rises[0]:
            index = 0
        elif pressure_rise >= pressure_rises[-1]:
            index = len(pressure_rises) - 2
        else:
            index = int(numpy.searchsorted(pressure_rises, pressure_rise) - 1)
        rise_low, rise_high = pressure_rises[index], pressure_rises[index + 1]
        flow_low, flow_high = mass_flows[index], mass_flows[index + 1]
        slope = (flow_high - flow_low) / (rise_high - rise_low) if rise_high != rise_low else 0.0
        return flow_low + slope * (pressure_rise - rise_low), slope
    return pump_curve


def parse_pump_curve_csv(path, density, ask_revolutions_per_minute=None):
    """Read a pump-curve CSV.  Returns {"rpm": [[pressure_rise_psid, flow_mL_per_min], ...]}.

    Header columns (case-insensitive):
      rpm                               (optional; ask_revolutions_per_minute() is called if absent)
      dp_psid | dp_psi | dp_bar | dp_kpa
      flow_ml_min | flow_l_min | flow_lpm | flow_gpm | flow_g_s
    """
    with open(path, newline='', encoding='utf-8-sig') as file_handle:
        rows = [row for row in csv.reader(file_handle) if row and any(cell.strip() for cell in row)]
    if len(rows) < 2:
        raise ValueError("CSV needs a header row and at least one data row.")
    header = [heading.strip().lower().replace(' ', '_').replace('/', '_') for heading in rows[0]]

    def find_column(*names):
        for name in names:
            if name in header:
                return header.index(name), name
        return None, None
    speed_column, _ = find_column('rpm', 'speed_rpm')
    pressure_column, pressure_name = find_column('dp_psid', 'dp_psi', 'dp_bar', 'dp_kpa', 'delta_p_psid', 'dp')
    flow_column, flow_name = find_column('flow_ml_min', 'flow_l_min', 'flow_lpm', 'flow_gpm', 'flow_g_s', 'flow')
    if pressure_column is None or flow_column is None:
        raise ValueError("Header must contain a pressure column (dp_psid, dp_bar, dp_kpa) "
                         "and a flow column (flow_ml_min, flow_lpm, flow_gpm, flow_g_s).")
    psid_per_unit = {'dp_psid': 1.0, 'dp_psi': 1.0, 'delta_p_psid': 1.0, 'dp': 1.0,
                     'dp_bar': 1e5 / PASCALS_PER_PSI, 'dp_kpa': 1e3 / PASCALS_PER_PSI}[pressure_name]
    milliliters_per_minute_per_unit = {
        'flow_ml_min': 1.0, 'flow': 1.0, 'flow_l_min': 1000.0, 'flow_lpm': 1000.0,
        'flow_gpm': CUBIC_METERS_PER_SECOND_PER_GALLON_PER_MINUTE / CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE,
        'flow_g_s': (KILOGRAMS_PER_GRAM / density) / CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE}[flow_name]
    fixed_speed = None
    if speed_column is None:
        fixed_speed = ask_revolutions_per_minute() if ask_revolutions_per_minute else None
        if fixed_speed is None:
            raise ValueError("No 'rpm' column and no RPM given.")
    curves = {}
    for row in rows[1:]:
        try:
            speed = float(row[speed_column]) if speed_column is not None else float(fixed_speed)
            pressure_rise = float(row[pressure_column]) * psid_per_unit
            flow = float(row[flow_column]) * milliliters_per_minute_per_unit
        except (ValueError, IndexError):
            continue
        curves.setdefault(str(int(round(speed))), []).append([pressure_rise, flow])
    for key in curves:
        curves[key].sort()
    if not curves:
        raise ValueError("No numeric data rows found.")
    return curves


# =============================================================================
# Network construction
# =============================================================================
class UnionFind:
    """Disjoint-set forest with path compression   [TARJAN]."""

    def __init__(self):
        self.parent = {}

    def find(self, item):
        self.parent.setdefault(item, item)
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:                       # path compression
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, first, second):
        first_root, second_root = self.find(first), self.find(second)
        if first_root != second_root:
            # keep fixed-boundary keys as roots for readability
            if second_root.startswith('__') or second_root.endswith(':reservoir'):
                first_root, second_root = second_root, first_root
            self.parent[second_root] = first_root


class NetworkElement:
    """A two-node network element.  Positive mass flow is first_node -> second_node.

    kind: 'resistance'  p_first - p_second = pressure_drop(mass_flow)
          'closed'      mass_flow = 0
          'fixed_flow'  mass_flow = fixed_mass_flow
          'pump_curve'  mass_flow = pump_curve(p_second - p_first)
    """
    __slots__ = ('kind', 'first_node', 'second_node', 'owner', 'pressure_drop_function',
                 'fixed_mass_flow', 'pump_curve', 'key')

    def __init__(self, kind, first_node, second_node, owner, key, pressure_drop_function=None,
                 fixed_mass_flow=0.0, pump_curve=None):
        self.kind, self.first_node, self.second_node = kind, first_node, second_node
        self.owner, self.key = owner, key
        self.pressure_drop_function = pressure_drop_function
        self.fixed_mass_flow, self.pump_curve = fixed_mass_flow, pump_curve


# Regularisation of |mass flow| in quadratic loss terms so that d(pressure drop)/d(mass flow)
# is never exactly zero at zero flow (keeps the Newton Jacobian non-singular)   [DENNIS]
MASS_FLOW_REGULARIZATION = 1e-7   # kg/s


def _signed_square(mass_flow):
    """mass_flow * |mass_flow|, smoothed at zero flow."""
    return mass_flow * math.sqrt(mass_flow * mass_flow + MASS_FLOW_REGULARIZATION ** 2)


def make_two_k_pressure_drop_function(turbulent_loss_coefficient, laminar_coefficient, diameter, density,
                                      viscosity, reverse_turbulent_loss_coefficient=None):
    """Pressure drop of a fitting referenced to `diameter`.

        K = K1 / Re + K_infinity                                    [HOOPER], cf. [BLEVINS]
        dP = K * rho * V^2 / 2                                      [CRANE] (resistance coefficient)
    With V = mdot / (rho * A) and Re = rho * V * D / mu:
        dP = K1 * mu * mdot / (2 rho A D)  +  K_infinity * mdot|mdot| / (2 rho A^2)
    An optional different K_infinity applies to reverse flow (e.g. entrance vs exit)."""
    area = math.pi * diameter * diameter / 4.0
    linear_coefficient = laminar_coefficient * viscosity / (2.0 * density * area * diameter)
    forward_quadratic_coefficient = turbulent_loss_coefficient / (2.0 * density * area * area)
    reverse_quadratic_coefficient = (forward_quadratic_coefficient if reverse_turbulent_loss_coefficient is None
                                     else reverse_turbulent_loss_coefficient / (2.0 * density * area * area))
    return lambda mass_flow: (linear_coefficient * mass_flow
                              + (forward_quadratic_coefficient if mass_flow >= 0 else reverse_quadratic_coefficient)
                              * _signed_square(mass_flow))


def make_orifice_pressure_drop_function(diameter, discharge_coefficient, density):
    """Equivalent sharp-edged orifice:
        mdot = Cd * A * sqrt(2 * rho * dP)                        [SUTTON], [CRANE]
    so  dP = mdot|mdot| / (2 * rho * Cd^2 * A^2)."""
    area = math.pi * diameter * diameter / 4.0
    quadratic_coefficient = 1.0 / (2.0 * density * discharge_coefficient ** 2 * area * area)
    return lambda mass_flow: quadratic_coefficient * _signed_square(mass_flow)


def area_change_quadratic_coefficient(upstream_area, downstream_area, density):
    """Coefficient c in dP = c * mdot|mdot| for flow from upstream_area into downstream_area.

    Sudden contraction:  K = 0.5 * (1 - beta^2), on the small (downstream) velocity   [CRANE]
    Sudden enlargement:  K = (1 - beta^2)^2,     on the small (upstream) velocity     [CRANE], [WHITE]
                         (Borda-Carnot loss), beta^2 = small area / large area."""
    if downstream_area < upstream_area:
        loss_coefficient = 0.5 * (1.0 - downstream_area / upstream_area)
        return loss_coefficient / (2.0 * density * downstream_area * downstream_area)
    loss_coefficient = (1.0 - upstream_area / downstream_area) ** 2
    return loss_coefficient / (2.0 * density * upstream_area * upstream_area)


def make_tube_pressure_drop_function(diameter, length, bend_count, bend_radius_ratio, roughness, density,
                                     viscosity, settings, end_area_change_coefficients=()):
    """Pressure drop of a straight tube with bends and optional end area changes.

    Friction (Darcy-Weisbach):  dP = f * (L/D) * rho * V^2 / 2              [CRANE], [WHITE]
        written as (f*Re) * mu * V * L / (2 D^2) with f from                [CHURCHILL]
    Bends (each 90°):           K = n(r/D) * f_T + K1/Re                    [CRANE], [HOOPER]
    End area changes:           see area_change_quadratic_coefficient()     [CRANE]"""
    area = math.pi * diameter * diameter / 4.0
    relative_roughness = roughness / diameter
    turbulent_friction = fully_turbulent_friction_factor(diameter, settings)
    bend_turbulent_coefficient = bend_count * crane_bend_multiplier(bend_radius_ratio) * turbulent_friction
    bend_laminar_coefficient = bend_count * settings.get('bend_laminar_coefficient', 800.0)
    bend_quadratic_coefficient = bend_turbulent_coefficient / (2.0 * density * area * area)
    bend_linear_coefficient = bend_laminar_coefficient * viscosity / (2.0 * density * area * diameter)
    end_coefficients = list(end_area_change_coefficients)

    def pressure_drop(mass_flow):
        reynolds_number = abs(mass_flow) * diameter / (area * viscosity)          # Re = rho V D / mu [WHITE]
        velocity = mass_flow / (density * area)
        drop = (churchill_friction_factor_times_reynolds(reynolds_number, relative_roughness)
                * viscosity * velocity * length / (2.0 * diameter * diameter))   # [CRANE] Darcy, [CHURCHILL]
        drop += bend_linear_coefficient * mass_flow + bend_quadratic_coefficient * _signed_square(mass_flow)
        for forward_coefficient, reverse_coefficient in end_coefficients:
            drop += (forward_coefficient if mass_flow >= 0 else reverse_coefficient) * _signed_square(mass_flow)
        return drop
    return pressure_drop


class Network:
    def __init__(self):
        self.union_find = UnionFind()
        self.elements = []
        self.fixed_pressures = {}   # raw node key -> pressure [Pa]
        self.issues = []            # (owner identifier, text)
        self.line_sections = {}     # line identifier -> list of (start fraction, end fraction,
        #                             element key, start node key, end node key, diameter)
        self.density = DEFAULT_FLUID['density']
        self.viscosity = DEFAULT_FLUID['viscosity_millipascal_seconds'] * 1e-3
        self.settings = dict(DEFAULT_SETTINGS)
        self.backfill = {}          # tank identifier -> (ullage pressure psia, COPV identifier)

    def node(self, key):
        return self.union_find.find(key)


def blocks_by_identifier(model):
    return {block['identifier']: block for block in model.get('blocks', [])}


def build_network(model):
    """Turn a schematic model into nodes and two-node elements."""
    fluid = model.get('fluid', DEFAULT_FLUID)
    settings = dict(DEFAULT_SETTINGS)
    settings.update(model.get('settings', {}))
    density = float(fluid['density'])
    viscosity = float(fluid['viscosity_millipascal_seconds']) * 1e-3
    network = Network()
    network.density, network.viscosity, network.settings = density, viscosity, settings
    network.backfill = ullage_backfill(model)
    blocks = blocks_by_identifier(model)
    lines = {line['identifier']: line for line in model.get('lines', [])}
    pending_elements = []   # (kind, first node, second node, owner, key, payload)

    # ---- line end points and taps -----------------------------------------------
    taps = {line_identifier: [] for line_identifier in lines}
    for line in lines.values():
        line_identifier = line['identifier']
        for end_name in ('start', 'end'):
            attachment = line.get(end_name)
            end_key = f"{line_identifier}:{end_name}"
            network.union_find.find(end_key)
            if not attachment:
                continue
            if attachment.get('kind') == 'port' and attachment.get('block') in blocks:
                network.union_find.union(f"{attachment['block']}:{attachment['port']}", end_key)
            elif (attachment.get('kind') == 'line' and attachment.get('line') in lines
                  and attachment['line'] != line_identifier):
                path_fraction = float(attachment.get('path_fraction', 0.0))
                if path_fraction <= 1e-6:
                    network.union_find.union(f"{attachment['line']}:start", end_key)
                elif path_fraction >= 1 - 1e-6:
                    network.union_find.union(f"{attachment['line']}:end", end_key)
                else:
                    taps[attachment['line']].append((path_fraction, end_key))

    # ---- line sections (a line is split where other lines tap into it) ------------
    for line in lines.values():
        line_identifier = line['identifier']
        properties = line.get('properties', {})
        diameter = max(float(properties.get('inner_diameter', 0.18)), 1e-4) * METERS_PER_INCH
        total_length = max(float(properties.get('length', 1.0)), 1e-3) * METERS_PER_INCH
        bend_count = max(float(properties.get('bend_count', 0)), 0.0)
        bend_radius_ratio = float(properties.get('bend_radius_ratio', 3.0))
        roughness = max(float(properties.get('roughness', 0.0)), 0.0) * METERS_PER_INCH
        line_area = math.pi * diameter * diameter / 4
        # area change where the tube meets a port of a different inner diameter   [CRANE]
        end_area_change = {'start': None, 'end': None}
        for end_name in ('start', 'end'):
            attachment = line.get(end_name)
            if attachment and attachment.get('kind') == 'port' and attachment.get('block') in blocks:
                block = blocks[attachment['block']]
                tubing_key = PORT_TUBING_PROPERTY.get(block['type'], {}).get(attachment['port'])
                if tubing_key:
                    port_diameter = tube_inner_diameter_inches(
                        block['properties'].get(tubing_key, DEFAULT_TUBE)) * METERS_PER_INCH
                    if abs(port_diameter - diameter) / diameter > settings['size_mismatch_tolerance']:
                        port_area = math.pi * port_diameter * port_diameter / 4
                        port_to_line = area_change_quadratic_coefficient(port_area, line_area, density)
                        line_to_port = area_change_quadratic_coefficient(line_area, port_area, density)
                        # positive line flow is start -> end
                        end_area_change[end_name] = ((port_to_line, line_to_port) if end_name == 'start'
                                                     else (line_to_port, port_to_line))
        sorted_taps = sorted(taps[line_identifier])
        cut_fractions = [0.0]
        cut_nodes = [f"{line_identifier}:start"]
        for path_fraction, end_key in sorted_taps:
            if path_fraction - cut_fractions[-1] < 1e-6:
                network.union_find.union(cut_nodes[-1], end_key)
                continue
            tap_node = f"{line_identifier}:tap{len(cut_fractions)}"
            cut_fractions.append(path_fraction)
            cut_nodes.append(tap_node)
            network.union_find.union(tap_node, end_key)
        cut_fractions.append(1.0)
        cut_nodes.append(f"{line_identifier}:end")
        sections = []
        section_count = len(cut_fractions) - 1
        for section_index in range(section_count):
            start_fraction, end_fraction = cut_fractions[section_index], cut_fractions[section_index + 1]
            fraction = end_fraction - start_fraction
            end_coefficients = []
            if section_index == 0 and end_area_change['start']:
                end_coefficients.append(end_area_change['start'])
            if section_index == section_count - 1 and end_area_change['end']:
                end_coefficients.append(end_area_change['end'])
            element_key = f"{line_identifier}#{section_index}"
            pressure_drop_function = make_tube_pressure_drop_function(
                diameter, total_length * fraction, bend_count * fraction, bend_radius_ratio, roughness,
                density, viscosity, settings, end_coefficients)
            pending_elements.append(('resistance', cut_nodes[section_index], cut_nodes[section_index + 1],
                                     ('line', line_identifier, section_index), element_key,
                                     pressure_drop_function))
            sections.append((start_fraction, end_fraction, element_key, cut_nodes[section_index],
                             cut_nodes[section_index + 1], diameter))
        network.line_sections[line_identifier] = sections

    # ---- blocks ----------------------------------------------------------------------
    orifice_discharge_coefficient = settings['orifice_discharge_coefficient']
    for block in blocks.values():
        identifier, block_type = block['identifier'], block['type']
        properties, state = block.get('properties', {}), normalized_state(block['type'], block.get('state', 0))
        inlet_node, outlet_node = f"{identifier}:inlet", f"{identifier}:outlet"
        for port_name in BLOCK_TYPES[block_type]['ports']:
            network.union_find.find(f"{identifier}:{port_name}")

        def add_orifice(first_node, second_node, diameter_inches, discharge_coefficient, key='main'):
            """Equivalent sharp-edged orifice element, or a closed element for zero diameter."""
            if diameter_inches and diameter_inches > 0 and discharge_coefficient > 0:
                pending_elements.append(('resistance', first_node, second_node, ('block', identifier, key),
                                         f"{identifier}#{key}",
                                         make_orifice_pressure_drop_function(diameter_inches * METERS_PER_INCH,
                                                                             discharge_coefficient, density)))
            else:
                pending_elements.append(('closed', first_node, second_node, ('block', identifier, key),
                                         f"{identifier}#{key}", None))

        if block_type == 'tank':
            reservoir_node = f"{identifier}:reservoir"
            # ullage pressure: the backfilling COPV's pressure if one is connected, else the tank state
            tank_pressure = (network.backfill[identifier][0] if identifier in network.backfill
                             else tank_pressure_psia(properties, state))
            network.fixed_pressures[reservoir_node] = tank_pressure * PASCALS_PER_PSI
            outlet_diameter = tube_inner_diameter_inches(
                properties.get('outlet_tubing_size', DEFAULT_TUBE)) * METERS_PER_INCH
            # bottom outlet: entrance loss leaving the tank, exit loss for back-flow   [CRANE]
            pending_elements.append(('resistance', reservoir_node, outlet_node, ('block', identifier, 'main'),
                                     f"{identifier}#main",
                                     make_two_k_pressure_drop_function(
                                         settings['tank_entrance_loss_coefficient'], 0.0, outlet_diameter,
                                         density, viscosity,
                                         reverse_turbulent_loss_coefficient=settings['tank_exit_loss_coefficient'])))
            # top inlet: positive flow is into the tank (exit loss), reverse is entrance loss   [CRANE]
            inlet_diameter = tube_inner_diameter_inches(
                properties.get('inlet_tubing_size', DEFAULT_TUBE)) * METERS_PER_INCH
            pending_elements.append(('resistance', inlet_node, reservoir_node, ('block', identifier, 'inlet'),
                                     f"{identifier}#inlet",
                                     make_two_k_pressure_drop_function(
                                         settings['tank_exit_loss_coefficient'], 0.0, inlet_diameter,
                                         density, viscosity,
                                         reverse_turbulent_loss_coefficient=settings['tank_entrance_loss_coefficient'])))
        elif block_type == 'latch_valve':
            diameter_key = 'open_orifice_diameter' if state == 0 else 'closed_orifice_diameter'
            add_orifice(inlet_node, outlet_node, float(properties.get(diameter_key, 0)),
                        float(properties.get('discharge_coefficient', 0.6)))
        elif block_type == 'filter':
            diameter_key = 'clean_orifice_diameter' if state == 0 else 'clogged_orifice_diameter'
            add_orifice(inlet_node, outlet_node, float(properties.get(diameter_key, 0)),
                        float(properties.get('discharge_coefficient', 0.6)))
        elif block_type == 'flowmeter':
            if state == 0:
                add_orifice(inlet_node, outlet_node, float(properties.get('orifice_diameter', 0)),
                            float(properties.get('discharge_coefficient', 0.6)))
            else:
                add_orifice(inlet_node, outlet_node, 0.0, 1.0)
        elif block_type == 'transducer':
            pass  # dead-end pressure tap: reads the pressure of the node it is connected to
        elif block_type in ('copv', 'fill_drain_valve'):
            pass  # COPV: sets tank ullage pressure (see ullage_backfill); fill & drain valve: dead ends
        elif block_type == 'fluid_coupling':
            if state == 0:
                add_orifice(inlet_node, outlet_node, float(properties.get('orifice_diameter', 0)),
                            float(properties.get('discharge_coefficient', 0.6)))
            else:
                add_orifice(inlet_node, outlet_node, 0.0, 1.0)
        elif block_type == 'pump':
            # ON: flow = curve(pressure rise) at the set speed, so the outlet pressure is where the
            # pump curve meets the system.  OFF, or ON with no curves: stopped-pump orifice.
            curve_function = None
            if state == 0:
                curve_function, curve_note = make_pump_curve_function_for_speed(
                    block.get('pump_curves', {}), float(properties.get('revolutions_per_minute', 0)), density)
                if curve_note:
                    network.issues.append((identifier, curve_note))
                if curve_function is None:
                    network.issues.append((identifier, "Pump is ON but has no curves (right-click → Pump "
                                                       "curves… → Import CSV); solved as stopped (orifice)."))
            if curve_function is not None:
                pending_elements.append(('pump_curve', inlet_node, outlet_node,
                                         ('block', identifier, 'main'), f"{identifier}#main", curve_function))
            else:
                add_orifice(inlet_node, outlet_node, float(properties.get('stopped_orifice_diameter', 0)),
                            float(properties.get('stopped_discharge_coefficient', 0.6)))
        elif block_type == 'thruster':
            if state == 0:
                # thruster = equivalent orifice discharging to 0 psia (hard coded)
                add_orifice(inlet_node, '__ambient__', float(properties.get('orifice_diameter', 0)),
                            float(properties.get('discharge_coefficient', 0.6)))
                network.fixed_pressures['__ambient__'] = 0.0
        elif block_type == 'orifice':
            if state == 0:
                add_orifice(inlet_node, outlet_node, float(properties.get('orifice_diameter', 0)),
                            orifice_discharge_coefficient)
            else:
                add_orifice(inlet_node, outlet_node, 0.0, 1.0)
        elif block_type == 'cross':
            # Four arms meeting at a centre node.  Crane TP-410 gives no separate cross value,
            # so the cross uses the tee coefficients: flow-through-run 20 f_T split over the
            # inlet and outlet arms, flow-through-branch 60 f_T for inlet arm + either branch
            # arm   [CRANE]; laminar terms split the same way   [HOOPER].
            center_node = f"{identifier}:center"
            run_multiplier = settings['tee_run_friction_multiplier']
            branch_multiplier = settings['tee_branch_friction_multiplier']
            run_laminar = settings['tee_run_laminar_coefficient']
            branch_laminar = settings['tee_branch_laminar_coefficient']
            arms = [('inlet', inlet_node, center_node, 'inlet_tubing_size', run_multiplier / 2, run_laminar / 2),
                    ('outlet', center_node, outlet_node, 'outlet_tubing_size', run_multiplier / 2, run_laminar / 2)]
            for branch_name in ('upper_branch', 'lower_branch'):
                arms.append((branch_name, center_node, f"{identifier}:{branch_name}", f"{branch_name}_tubing_size",
                             branch_multiplier - run_multiplier / 2, branch_laminar - run_laminar / 2))
            for arm_name, first_node, second_node, tubing_key, friction_multiplier, laminar_coefficient in arms:
                arm_diameter = tube_inner_diameter_inches(properties.get(tubing_key, DEFAULT_TUBE)) * METERS_PER_INCH
                # cross arm: K_turbulent = n * f_T                                           [CRANE]
                # plus laminar term K_1 / Re (two-K method)                              [HOOPER]
                pending_elements.append(('resistance', first_node, second_node, ('block', identifier, arm_name),
                                         f"{identifier}#{arm_name}",
                                         make_two_k_pressure_drop_function(
                                             friction_multiplier * fully_turbulent_friction_factor(arm_diameter, settings),
                                             laminar_coefficient, arm_diameter, density, viscosity)))
        elif block_type == 'tee':
            # Three arms meeting at a centre node.  Crane flow-through-run (20 f_T) and
            # flow-through-branch (60 f_T) are split so inlet+outlet arms = run and
            # inlet+branch arms = branch   [CRANE]; laminar terms likewise   [HOOPER].
            center_node = f"{identifier}:center"
            run_multiplier = settings['tee_run_friction_multiplier']
            branch_multiplier = settings['tee_branch_friction_multiplier']
            run_laminar = settings['tee_run_laminar_coefficient']
            branch_laminar = settings['tee_branch_laminar_coefficient']
            arms = [('inlet', inlet_node, center_node, 'inlet_tubing_size', run_multiplier / 2, run_laminar / 2),
                    ('outlet', center_node, outlet_node, 'outlet_tubing_size', run_multiplier / 2, run_laminar / 2),
                    ('branch', center_node, f"{identifier}:branch", 'branch_tubing_size',
                     branch_multiplier - run_multiplier / 2, branch_laminar - run_laminar / 2)]
            for arm_name, first_node, second_node, tubing_key, friction_multiplier, laminar_coefficient in arms:
                arm_diameter = tube_inner_diameter_inches(properties.get(tubing_key, DEFAULT_TUBE)) * METERS_PER_INCH
                # tee arm: K_turbulent = n * f_T (run or branch multiplier)                   [CRANE]
                # plus laminar term K_1 / Re (two-K method)                              [HOOPER]
                pending_elements.append(('resistance', first_node, second_node, ('block', identifier, arm_name),
                                         f"{identifier}#{arm_name}",
                                         make_two_k_pressure_drop_function(
                                             friction_multiplier * fully_turbulent_friction_factor(arm_diameter, settings),
                                             laminar_coefficient, arm_diameter, density, viscosity)))
        elif block_type == 'elbow':
            inlet_diameter = tube_inner_diameter_inches(
                properties.get('inlet_tubing_size', DEFAULT_TUBE)) * METERS_PER_INCH
            outlet_diameter = tube_inner_diameter_inches(
                properties.get('outlet_tubing_size', DEFAULT_TUBE)) * METERS_PER_INCH
            # standard 90° elbow K = 30 f_T [CRANE] with laminar term K1/Re [HOOPER]
            elbow_function = make_two_k_pressure_drop_function(
                settings['elbow_friction_multiplier'] * fully_turbulent_friction_factor(inlet_diameter, settings),
                settings['elbow_laminar_coefficient'], inlet_diameter, density, viscosity)
            if abs(outlet_diameter - inlet_diameter) / inlet_diameter > settings['size_mismatch_tolerance']:
                inlet_area = math.pi * inlet_diameter * inlet_diameter / 4
                outlet_area = math.pi * outlet_diameter * outlet_diameter / 4
                forward = area_change_quadratic_coefficient(inlet_area, outlet_area, density)
                reverse = area_change_quadratic_coefficient(outlet_area, inlet_area, density)

                def make_combined(base_function, forward_coefficient, reverse_coefficient):
                    return lambda mass_flow: (base_function(mass_flow)
                                              + (forward_coefficient if mass_flow >= 0 else reverse_coefficient)
                                              * _signed_square(mass_flow))
                pressure_drop_function = make_combined(elbow_function, forward, reverse)
            else:
                pressure_drop_function = elbow_function
            pending_elements.append(('resistance', inlet_node, outlet_node, ('block', identifier, 'main'),
                                     f"{identifier}#main", pressure_drop_function))

    # finalise elements with merged node keys
    for kind, first_node, second_node, owner, key, payload in pending_elements:
        element = NetworkElement(kind, network.union_find.find(first_node), network.union_find.find(second_node),
                                 owner, key)
        if kind == 'resistance':
            element.pressure_drop_function = payload
        elif kind == 'fixed_flow':
            element.fixed_mass_flow = payload
        elif kind == 'pump_curve':
            element.pump_curve = payload
        network.elements.append(element)
    network.fixed_pressures = {network.union_find.find(key): value
                               for key, value in network.fixed_pressures.items()}
    return network


# =============================================================================
# Solver
# =============================================================================
PRESSURE_SCALE = 1e5      # Pa      - unknowns and residuals are scaled to order one
MASS_FLOW_SCALE = 1e-3    # kg/s


def _central_difference_derivative(function, value):
    """Central finite difference f'(x) ~ (f(x+h) - f(x-h)) / 2h   [DENNIS]."""
    step = max(1e-9, 1e-6 * abs(value))
    return (function(value + step) - function(value - step)) / (2 * step)


def solve_network(network, initial_guess=None, tolerance=1e-10, maximum_iterations=150):
    """Solve node pressures and element mass flows simultaneously   [TODINI].

    Unknowns: pressure of every free node and mass flow of every active element.
    Equations: mass conservation at every free node, plus each element's law.
    Method: Newton-Raphson with a backtracking line search   [DENNIS], started
    from a linearised pre-solve.  Nodes in regions with no pressure boundary get NaN.

    Returns dict(converged, message, iterations, pressures={node: Pa},
                 mass_flows={element key: kg/s}, issues=[...]).
    """
    elements = network.elements
    nodes = set()
    for element in elements:
        nodes.add(element.first_node)
        nodes.add(element.second_node)
    for key in list(network.union_find.parent.keys()):
        nodes.add(network.union_find.find(key))
    nodes |= set(network.fixed_pressures)

    issues = []
    # pressure-transmitting connectivity: which regions touch a fixed-pressure node?
    for _ in range(3):
        connectivity = UnionFind()
        for node in nodes:
            connectivity.find(node)
        for element in elements:
            if element.kind in ('resistance', 'pump_curve'):
                connectivity.union(element.first_node, element.second_node)
        anchored_roots = {connectivity.find(node) for node in network.fixed_pressures}
        changed = False
        for element in elements:
            if element.kind == 'fixed_flow' and abs(element.fixed_mass_flow) > 0:
                first_anchored = connectivity.find(element.first_node) in anchored_roots
                second_anchored = connectivity.find(element.second_node) in anchored_roots
                if not (first_anchored and second_anchored):
                    side = "inlet" if not first_anchored else "outlet"
                    issues.append((element.owner[1], f"Pump {side} has no open path to a pressure "
                                                     f"boundary (dead-headed/starved). Flow forced to 0."))
                    element.kind = 'closed'
                    changed = True
        if not changed:
            break
    floating_nodes = {node for node in nodes if connectivity.find(node) not in anchored_roots}

    free_nodes = sorted(node for node in nodes
                        if node not in network.fixed_pressures and node not in floating_nodes)
    node_index = {node: index for index, node in enumerate(free_nodes)}
    active_elements = [element for element in elements if element.kind != 'closed'
                       and element.first_node not in floating_nodes and element.second_node not in floating_nodes]
    node_count, element_count = len(free_nodes), len(active_elements)
    unknown_count = node_count + element_count
    result = {'converged': True, 'message': '', 'iterations': 0, 'pressures': {}, 'mass_flows': {},
              'issues': issues}
    for node, value in network.fixed_pressures.items():
        result['pressures'][node] = value
    for node in floating_nodes:
        result['pressures'][node] = float('nan')
    for element in elements:
        result['mass_flows'][element.key] = 0.0
    if unknown_count == 0:
        result['message'] = 'Nothing to solve.'
        return result

    def pressure_of(unknowns, node):
        if node in node_index:
            return unknowns[node_index[node]] * PRESSURE_SCALE
        return network.fixed_pressures.get(node, 0.0)

    def residuals_and_jacobian(unknowns, linearised=False, linear_resistances=None):
        residuals = numpy.zeros(unknown_count)
        jacobian = numpy.zeros((unknown_count, unknown_count))
        for element_index, element in enumerate(active_elements):
            row = node_count + element_index
            column = node_count + element_index
            mass_flow = unknowns[column] * MASS_FLOW_SCALE
            first_pressure = pressure_of(unknowns, element.first_node)
            second_pressure = pressure_of(unknowns, element.second_node)
            first_index = node_index.get(element.first_node)
            second_index = node_index.get(element.second_node)
            if element.kind == 'resistance':
                # element law: p_first - p_second - dP(mdot) = 0
                if linearised:
                    drop, slope = linear_resistances[element_index] * mass_flow, linear_resistances[element_index]
                else:
                    drop = element.pressure_drop_function(mass_flow)
                    slope = _central_difference_derivative(element.pressure_drop_function, mass_flow)
                residuals[row] = (first_pressure - second_pressure - drop) / PRESSURE_SCALE
                if first_index is not None:
                    jacobian[row, first_index] += 1.0
                if second_index is not None:
                    jacobian[row, second_index] -= 1.0
                jacobian[row, column] = -slope * MASS_FLOW_SCALE / PRESSURE_SCALE
            elif element.kind == 'fixed_flow':
                # element law: mdot - mdot_set = 0
                residuals[row] = (mass_flow - element.fixed_mass_flow) / MASS_FLOW_SCALE
                jacobian[row, column] = 1.0
            elif element.kind == 'pump_curve':
                # element law: mdot - curve(p_second - p_first) = 0
                if linearised:
                    curve_flow, slope = element.pump_curve(0.0)
                    curve_flow = curve_flow + slope * (second_pressure - first_pressure)
                else:
                    curve_flow, slope = element.pump_curve(second_pressure - first_pressure)
                residuals[row] = (mass_flow - curve_flow) / MASS_FLOW_SCALE
                jacobian[row, column] = 1.0
                if first_index is not None:
                    jacobian[row, first_index] += slope * PRESSURE_SCALE / MASS_FLOW_SCALE
                if second_index is not None:
                    jacobian[row, second_index] -= slope * PRESSURE_SCALE / MASS_FLOW_SCALE
            # mass conservation at the element's end nodes: sum(in) - sum(out) = 0
            if first_index is not None:
                residuals[first_index] -= unknowns[column]
                jacobian[first_index, column] -= 1.0
            if second_index is not None:
                residuals[second_index] += unknowns[column]
                jacobian[second_index, column] += 1.0
        return residuals, jacobian

    # linearised pre-solve (each resistance replaced by dP = R * mdot) for a robust start
    unknowns = numpy.zeros(unknown_count)
    fixed_values = list(network.fixed_pressures.values())
    starting_pressure = (max(fixed_values) if fixed_values else 0.0) * 0.5
    unknowns[:node_count] = starting_pressure / PRESSURE_SCALE
    reference_mass_flow = 1e-3
    linear_resistances = [max(element.pressure_drop_function(reference_mass_flow) / reference_mass_flow, 1.0)
                          if element.kind == 'resistance' else 0.0 for element in active_elements]
    try:
        residuals, jacobian = residuals_and_jacobian(unknowns, linearised=True,
                                                     linear_resistances=linear_resistances)
        # one Newton step on the linearised system: x <- x - J^-1 F(x)                     [DENNIS]
        unknowns = unknowns - numpy.linalg.lstsq(jacobian, residuals, rcond=None)[0]
    except numpy.linalg.LinAlgError:
        pass
    if initial_guess:
        for node, index in node_index.items():
            value = initial_guess.get('pressures', {}).get(node)
            if value is not None and not math.isnan(value):
                unknowns[index] = value / PRESSURE_SCALE
        for element_index, element in enumerate(active_elements):
            value = initial_guess.get('mass_flows', {}).get(element.key)
            if value is not None:
                unknowns[node_count + element_index] = value / MASS_FLOW_SCALE

    # Newton-Raphson with backtracking line search   [DENNIS]
    residuals, jacobian = residuals_and_jacobian(unknowns)
    residual_norm = float(numpy.linalg.norm(residuals))
    iteration = 0
    converged = False
    for iteration in range(1, maximum_iterations + 1):
        if numpy.max(numpy.abs(residuals)) < tolerance:
            converged = True
            break
        try:
            step = numpy.linalg.solve(jacobian, -residuals)               # J * step = -F
        except numpy.linalg.LinAlgError:
            step = numpy.linalg.lstsq(jacobian, -residuals, rcond=None)[0]
        step_fraction = 1.0
        # backtracking line search with sufficient-decrease (Armijo-type) test:
        #   ||F(x + t s)|| < (1 - 1e-4 t) ||F(x)||, halving t each try                   [DENNIS]
        while True:
            trial_unknowns = unknowns + step_fraction * step
            trial_residuals, trial_jacobian = residuals_and_jacobian(trial_unknowns)
            trial_norm = float(numpy.linalg.norm(trial_residuals))
            if trial_norm < residual_norm * (1 - 1e-4 * step_fraction) or step_fraction < 1e-6:
                break
            step_fraction *= 0.5
        unknowns, residuals, jacobian, residual_norm = (trial_unknowns, trial_residuals,
                                                        trial_jacobian, trial_norm)
    else:
        converged = numpy.max(numpy.abs(residuals)) < 1e-6

    result['converged'] = bool(converged)
    result['iterations'] = iteration
    result['residual'] = float(numpy.max(numpy.abs(residuals))) if unknown_count else 0.0
    result['message'] = (f"Converged in {iteration} iterations" if converged
                         else f"Did not converge (residual {result['residual']:.2e})")
    for node, index in node_index.items():
        result['pressures'][node] = float(unknowns[index] * PRESSURE_SCALE)
    for element_index, element in enumerate(active_elements):
        result['mass_flows'][element.key] = float(unknowns[node_count + element_index] * MASS_FLOW_SCALE)
    return result


# =============================================================================
# Model-level solve and post-processing
# =============================================================================
def thruster_nominal_mass_flow(properties, density):
    """Flow at nominal inlet pressure into a 0 psia exhaust:
        mdot = Cd * A * sqrt(2 * rho * dP)                        [SUTTON], [CRANE]"""
    diameter = float(properties.get('orifice_diameter', 0)) * METERS_PER_INCH
    area = math.pi * diameter * diameter / 4
    pressure_drop = max(float(properties.get('nominal_inlet_pressure', 0)), 0) * PASCALS_PER_PSI
    return float(properties.get('discharge_coefficient', 0.6)) * area * math.sqrt(2 * density * pressure_drop)


def psia_from_pascals(pressure_pascals):
    return pressure_pascals / PASCALS_PER_PSI


def solve_model(model, initial_guess=None):
    """Solve a schematic model.  Returns a results dictionary for the GUI."""
    network = build_network(model)
    solution = solve_network(network, initial_guess=initial_guess)
    density = network.density

    def pressure_at(key):
        return solution['pressures'].get(network.union_find.find(key), float('nan'))

    def mass_flow_of(key):
        return solution['mass_flows'].get(key, 0.0)

    results = {'converged': solution['converged'], 'message': solution['message'],
               'iterations': solution['iterations'], 'blocks': {}, 'lines': {}, 'warnings': [],
               'initial_guess': {'pressures': solution['pressures'], 'mass_flows': solution['mass_flows']},
               'density': density}
    results['warnings'].extend(network.issues)
    results['warnings'].extend(solution['issues'])

    issues_by_identifier = {}
    for identifier, text in network.issues:
        issues_by_identifier.setdefault(identifier, []).append(text)
    blocks_warned_for_negative_pressure = set()
    for block in model.get('blocks', []):
        identifier, block_type = block['identifier'], block['type']
        properties, state = block.get('properties', {}), normalized_state(block['type'], block.get('state', 0))
        block_result = {'ports': {}, 'calculated': {}, 'label': ''}
        for port_name in BLOCK_TYPES[block_type]['ports']:
            block_result['ports'][port_name] = pressure_at(f"{identifier}:{port_name}")
        warnings = []
        calculated = block_result['calculated']
        inlet_pressure = block_result['ports'].get('inlet', float('nan'))
        outlet_pressure = block_result['ports'].get('outlet', float('nan'))
        mass_flow = mass_flow_of(f"{identifier}#main")
        if block_type in ('latch_valve', 'fluid_coupling', 'filter', 'flowmeter', 'orifice', 'elbow'):
            calculated['calculated_inlet_pressure'] = psia_from_pascals(inlet_pressure)
            calculated['calculated_outlet_pressure'] = psia_from_pascals(outlet_pressure)
            calculated['calculated_pressure_drop'] = psia_from_pascals(inlet_pressure - outlet_pressure)
            calculated['calculated_mass_flow'] = mass_flow / KILOGRAMS_PER_GRAM
            block_result['mass_flow'] = mass_flow
            block_result['pressure_drop'] = inlet_pressure - outlet_pressure
            block_result['label'] = (f"{format_number(mass_flow / KILOGRAMS_PER_GRAM, 3)} g/s · ΔP "
                                     f"{format_number(psia_from_pascals(inlet_pressure - outlet_pressure), 3)} psid")
            if block_type == 'orifice':
                discharge_coefficient = model.get('settings', {}).get(
                    'orifice_discharge_coefficient', DEFAULT_SETTINGS['orifice_discharge_coefficient'])
                calculated['calculated_discharge_coefficient'] = discharge_coefficient
                calculated['calculated_flow_coefficient'] = flow_coefficient_from_orifice_diameter(
                    float(properties.get('orifice_diameter', 0)), discharge_coefficient)
            maximum_pressure = properties.get('maximum_pressure')
            if maximum_pressure is not None:
                highest = max(value for value in (inlet_pressure, outlet_pressure, -1e30) if not math.isnan(value))
                if highest > float(maximum_pressure) * PASCALS_PER_PSI:
                    warnings.append(f"Pressure {format_number(psia_from_pascals(highest))} psia exceeds maximum "
                                    f"{format_number(float(maximum_pressure))} psia")
        elif block_type == 'tank':
            backfill_pressure, backfill_copv = network.backfill.get(identifier, (None, None))
            tank_pressure = backfill_pressure if backfill_copv else tank_pressure_psia(properties, state)
            if backfill_copv:
                block_result['backfill_copv'] = backfill_copv
                if tank_pressure > float(properties.get('maximum_pressure', 1e9)):
                    warnings.append(f"COPV ullage pressure {format_number(tank_pressure)} psia exceeds the tank "
                                    f"maximum {format_number(float(properties.get('maximum_pressure', 0)))} psia")
            inflow = mass_flow_of(f"{identifier}#inlet")
            net_outflow = mass_flow - inflow
            calculated['calculated_pressure'] = tank_pressure
            calculated['calculated_inflow'] = inflow / KILOGRAMS_PER_GRAM
            calculated['calculated_net_outflow'] = net_outflow / KILOGRAMS_PER_GRAM
            propellant_mass = float(properties.get('volume', 0)) * CUBIC_METERS_PER_LITER * density
            # time to empty = contained mass / net outflow (steady state)
            calculated['calculated_time_to_empty'] = (propellant_mass / net_outflow if net_outflow > 1e-12
                                                      else float('nan'))
            block_result['mass_flow'] = mass_flow
            block_result['label'] = (f"{format_number(tank_pressure, 4)} psia{' (COPV)' if backfill_copv else ''} · out "
                                     f"{format_number(mass_flow / KILOGRAMS_PER_GRAM, 3)} g/s")
            if abs(inflow) > 1e-9:
                block_result['label'] += f" · in {format_number(inflow / KILOGRAMS_PER_GRAM, 3)} g/s"
        elif block_type == 'copv':
            delivered = copv_delivered_pressure(properties, state)
            served_tanks = [tank for tank, (_, source) in network.backfill.items() if source == identifier]
            backfill_volume_flow = sum(mass_flow_of(f"{tank}#main") - mass_flow_of(f"{tank}#inlet")
                                       for tank in served_tanks) / density
            calculated['calculated_delivered_pressure'] = delivered if delivered is not None else float('nan')
            calculated['calculated_tanks_backfilled'] = len(served_tanks)
            # ullage gas must replace the propellant volume leaving the tank (incompressible liquid)
            calculated['calculated_backfill_volume_flow'] = (backfill_volume_flow
                                                             / CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE)
            if delivered is None:
                block_result['label'] = 'isolated'
            else:
                block_result['label'] = (f"ullage {format_number(delivered, 4)} psia → "
                                         f"{len(served_tanks)} tank{'s' if len(served_tanks) != 1 else ''}")
            bottle_pressure = float(properties.get('bottle_pressure', 0.0))
            setpoint = float(properties.get('ullage_pressure_setpoint', 0.0) or 0.0)
            if bottle_pressure > float(properties.get('maximum_pressure', 1e9)):
                warnings.append(f"Bottle pressure {format_number(bottle_pressure)} psia exceeds MEOP "
                                f"{format_number(float(properties.get('maximum_pressure', 0)))} psia")
            if delivered is not None and 0 < bottle_pressure < setpoint:
                warnings.append(f"Bottle pressure {format_number(bottle_pressure)} psia is below the setpoint "
                                f"{format_number(setpoint)} psia; delivering bottle pressure")
            if delivered is not None and not served_tanks:
                warnings.append("Not connected to a tank inlet (top port) through open lines and components")
        elif block_type == 'fill_drain_valve':
            calculated['calculated_inlet_pressure'] = psia_from_pascals(inlet_pressure)
            calculated['calculated_outlet_pressure'] = psia_from_pascals(outlet_pressure)
            block_result['label'] = ('dead end' if math.isnan(inlet_pressure)
                                     else f"dead end · {format_number(psia_from_pascals(inlet_pressure), 4)} psia")
            maximum_pressure = float(properties.get('maximum_pressure', 1e9))
            for port_pressure in (inlet_pressure, outlet_pressure):
                if not math.isnan(port_pressure) and port_pressure > maximum_pressure * PASCALS_PER_PSI:
                    warnings.append(f"Pressure {format_number(psia_from_pascals(port_pressure))} psia exceeds "
                                    f"maximum {format_number(maximum_pressure)} psia")
                    break
        elif block_type == 'transducer':
            calculated['calculated_pressure'] = psia_from_pascals(inlet_pressure)
            block_result['reading'] = psia_from_pascals(inlet_pressure) if state == 0 else float('nan')
            block_result['label'] = ''
        elif block_type == 'pump':
            calculated['calculated_inlet_pressure'] = psia_from_pascals(inlet_pressure)
            calculated['calculated_outlet_pressure'] = psia_from_pascals(outlet_pressure)
            calculated['calculated_pressure_rise'] = psia_from_pascals(outlet_pressure - inlet_pressure)
            calculated['calculated_mass_flow'] = mass_flow / KILOGRAMS_PER_GRAM
            calculated['calculated_volumetric_flow'] = (mass_flow / density
                                                        / CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE)
            block_result['mass_flow'] = mass_flow
            block_result['pressure_drop'] = inlet_pressure - outlet_pressure
            block_result['label'] = (f"{format_number(mass_flow / KILOGRAMS_PER_GRAM, 3)} g/s · rise "
                                     f"{format_number(psia_from_pascals(outlet_pressure - inlet_pressure), 4)} psid")
            if state == 0:
                limits = (
                    (inlet_pressure < float(properties.get('minimum_inlet_pressure', -1e9)) * PASCALS_PER_PSI,
                     f"Inlet {format_number(psia_from_pascals(inlet_pressure))} psia below minimum "
                     f"{format_number(float(properties.get('minimum_inlet_pressure', 0)))}"),
                    (inlet_pressure > float(properties.get('maximum_inlet_pressure', 1e9)) * PASCALS_PER_PSI,
                     f"Inlet {format_number(psia_from_pascals(inlet_pressure))} psia above maximum "
                     f"{format_number(float(properties.get('maximum_inlet_pressure', 0)))}"),
                    (outlet_pressure > float(properties.get('maximum_outlet_pressure', 1e9)) * PASCALS_PER_PSI,
                     f"Outlet {format_number(psia_from_pascals(outlet_pressure))} psia above maximum "
                     f"{format_number(float(properties.get('maximum_outlet_pressure', 0)))}"),
                    (outlet_pressure - inlet_pressure
                     > float(properties.get('maximum_differential_pressure', 1e9)) * PASCALS_PER_PSI,
                     f"ΔP {format_number(psia_from_pascals(outlet_pressure - inlet_pressure))} psid above maximum "
                     f"{format_number(float(properties.get('maximum_differential_pressure', 0)))}"))
                warnings.extend(text for violated, text in limits if violated)
        elif block_type == 'thruster':
            nominal_mass_flow = thruster_nominal_mass_flow(properties, density)
            specific_impulse = float(properties.get('specific_impulse', 0))
            calculated['calculated_outlet_pressure'] = 0.0
            calculated['calculated_inlet_pressure'] = psia_from_pascals(inlet_pressure)
            calculated['calculated_mass_flow'] = mass_flow / KILOGRAMS_PER_GRAM
            calculated['calculated_thrust'] = mass_flow * specific_impulse * STANDARD_GRAVITY   # F = mdot Isp g0 [SUTTON]
            calculated['calculated_nominal_mass_flow'] = nominal_mass_flow / KILOGRAMS_PER_GRAM
            calculated['calculated_nominal_thrust'] = nominal_mass_flow * specific_impulse * STANDARD_GRAVITY  # [SUTTON]
            block_result['mass_flow'] = mass_flow
            block_result['nominal_mass_flow'] = nominal_mass_flow
            block_result['label'] = (f"{format_number(mass_flow / KILOGRAMS_PER_GRAM, 3)} g/s · "
                                     f"{format_number(calculated['calculated_thrust'], 3)} N")
            if state == 0 and inlet_pressure > float(properties.get('maximum_inlet_pressure', 1e9)) * PASCALS_PER_PSI:
                warnings.append(f"Inlet {format_number(psia_from_pascals(inlet_pressure))} psia above maximum "
                                f"{format_number(float(properties.get('maximum_inlet_pressure', 0)))}")
        elif block_type == 'tee':
            calculated['calculated_center_pressure'] = psia_from_pascals(pressure_at(f"{identifier}:center"))
            calculated['calculated_inlet_flow'] = mass_flow_of(f"{identifier}#inlet") / KILOGRAMS_PER_GRAM
            calculated['calculated_outlet_flow'] = mass_flow_of(f"{identifier}#outlet") / KILOGRAMS_PER_GRAM
            calculated['calculated_branch_flow'] = mass_flow_of(f"{identifier}#branch") / KILOGRAMS_PER_GRAM
            block_result['label'] = (f"in {format_number(calculated['calculated_inlet_flow'], 3)} / out "
                                     f"{format_number(calculated['calculated_outlet_flow'], 3)} / branch "
                                     f"{format_number(calculated['calculated_branch_flow'], 3)} g/s")
        elif block_type == 'cross':
            calculated['calculated_center_pressure'] = psia_from_pascals(pressure_at(f"{identifier}:center"))
            calculated['calculated_inlet_flow'] = mass_flow_of(f"{identifier}#inlet") / KILOGRAMS_PER_GRAM
            calculated['calculated_outlet_flow'] = mass_flow_of(f"{identifier}#outlet") / KILOGRAMS_PER_GRAM
            calculated['calculated_upper_branch_flow'] = mass_flow_of(f"{identifier}#upper_branch") / KILOGRAMS_PER_GRAM
            calculated['calculated_lower_branch_flow'] = mass_flow_of(f"{identifier}#lower_branch") / KILOGRAMS_PER_GRAM
            block_result['label'] = (f"in {format_number(calculated['calculated_inlet_flow'], 3)} / out "
                                     f"{format_number(calculated['calculated_outlet_flow'], 3)} / up "
                                     f"{format_number(calculated['calculated_upper_branch_flow'], 3)} / down "
                                     f"{format_number(calculated['calculated_lower_branch_flow'], 3)} g/s")
        for port_name, port_pressure in block_result['ports'].items():
            if (not math.isnan(port_pressure) and port_pressure < -1.0
                    and identifier not in blocks_warned_for_negative_pressure):
                warnings.append(f"Negative absolute pressure at {port_name} "
                                f"({format_number(psia_from_pascals(port_pressure))} psia) — "
                                "flow demand cannot be met (cavitation).")
                blocks_warned_for_negative_pressure.add(identifier)
        for warning in warnings:
            results['warnings'].append((identifier, warning))
        # network issues (already in results['warnings']) are also shown on the block
        block_result['warnings'] = issues_by_identifier.get(identifier, []) + warnings
        results['blocks'][identifier] = block_result

    for line in model.get('lines', []):
        line_identifier = line['identifier']
        sections = []
        for (start_fraction, end_fraction, element_key, start_node, end_node,
             diameter) in network.line_sections.get(line_identifier, []):
            mass_flow = mass_flow_of(element_key)
            start_pressure, end_pressure = pressure_at(start_node), pressure_at(end_node)
            area = math.pi * diameter * diameter / 4                       # A = pi D^2 / 4
            # V = mdot / (rho A);  Re = rho V D / mu = |mdot| D / (A mu)                  [WHITE]
            sections.append({'start_fraction': start_fraction, 'end_fraction': end_fraction,
                             'mass_flow': mass_flow, 'start_pressure': start_pressure,
                             'end_pressure': end_pressure, 'pressure_drop': start_pressure - end_pressure,
                             'velocity': mass_flow / (density * area),
                             'reynolds_number': abs(mass_flow) * diameter / (area * network.viscosity)})
        line_result = {'sections': sections}
        if sections:
            first_section = sections[0]
            line_result['calculated'] = {
                'calculated_mass_flow': first_section['mass_flow'] / KILOGRAMS_PER_GRAM,
                'calculated_pressure_drop': (psia_from_pascals(first_section['start_pressure']
                                                               - sections[-1]['end_pressure'])
                                             * (1 if first_section['mass_flow'] >= 0 else -1)),
                'calculated_velocity': first_section['velocity'],
                'calculated_reynolds_number': first_section['reynolds_number']}
        results['lines'][line_identifier] = line_result
    return results


# =============================================================================
# Scenarios
# =============================================================================
def make_snapshot(model, name):
    return {'name': name,
            'blocks': {block['identifier']: {'state': int(block.get('state', 0)),
                                             'properties': copy.deepcopy(block.get('properties', {}))}
                       for block in model.get('blocks', [])},
            'lines': {line['identifier']: {'properties': copy.deepcopy(line.get('properties', {}))}
                      for line in model.get('lines', [])}}


def apply_snapshot(model, snapshot):
    """Apply a scenario snapshot to a model in place (unknown identifiers ignored)."""
    for block in model.get('blocks', []):
        saved = snapshot.get('blocks', {}).get(block['identifier'])
        if saved:
            block['state'] = normalized_state(block['type'], saved.get('state', 0))
            block['properties'].update(copy.deepcopy(saved.get('properties', {})))
    for line in model.get('lines', []):
        saved = snapshot.get('lines', {}).get(line['identifier'])
        if saved:
            line['properties'].update(copy.deepcopy(saved.get('properties', {})))
    return model


def find_block(model, identifier):
    for block in model.get('blocks', []):
        if block['identifier'] == identifier:
            return block
    return None


def thruster_demand(model, density):
    """Total nominal flow of all ON thrusters, and the per-thruster list."""
    total, per_thruster = 0.0, []
    for block in model.get('blocks', []):
        if block['type'] == 'thruster' and is_on(block):
            nominal = thruster_nominal_mass_flow(block['properties'], density)
            total += nominal
            per_thruster.append((block['identifier'], nominal))
    return total, per_thruster


def limit_violations(results, ignore=()):
    return [(identifier, text) for identifier, text in results['warnings'] if identifier not in ignore]


def evaluate_orifice(model, orifice_identifier, diameter_inches, initial_guess=None):
    """Solve the model with the given orifice diameter and summarise the flow balance
    pump flow = thruster flow + orifice back flow."""
    trial_model = copy.deepcopy(model)
    orifice_block = find_block(trial_model, orifice_identifier)
    orifice_block['properties']['orifice_diameter'] = float(diameter_inches)
    orifice_block['state'] = 0
    results = solve_model(trial_model, initial_guess=initial_guess)
    density = results['density']
    demand, _ = thruster_demand(trial_model, density)
    thruster_flow = sum(results['blocks'][block['identifier']].get('mass_flow', 0.0)
                        for block in trial_model['blocks'] if block['type'] == 'thruster' and is_on(block))
    pump_flow = sum(results['blocks'][block['identifier']].get('mass_flow', 0.0)
                    for block in trial_model['blocks'] if block['type'] == 'pump' and is_on(block))
    orifice_flow = results['blocks'][orifice_identifier].get('mass_flow', 0.0)
    relative_error = (thruster_flow / demand - 1.0) if demand > 0 else float('nan')
    pump_results = [results['blocks'][block['identifier']] for block in trial_model['blocks']
                    if block['type'] == 'pump' and is_on(block)]
    return {'diameter': diameter_inches, 'relative_error': relative_error, 'thruster_flow': thruster_flow,
            'demand': demand, 'pump_flow': pump_flow, 'orifice_flow': orifice_flow,
            'converged': results['converged'], 'violations': limit_violations(results), 'results': results,
            'pump_pressure_rise': max([pump['calculated']['calculated_pressure_rise'] for pump in pump_results],
                                      default=float('nan')),
            'pump_outlet_pressure': max([pump['calculated']['calculated_outlet_pressure'] for pump in pump_results],
                                        default=float('nan'))}


ORIFICE_DIAMETER_LOWER_BOUND, ORIFICE_DIAMETER_UPPER_BOUND = 1e-4, 0.5   # inches


def ideal_orifice_for_model(model, orifice_identifier, lower_bound=ORIFICE_DIAMETER_LOWER_BOUND,
                            upper_bound=ORIFICE_DIAMETER_UPPER_BOUND):
    """Diameter making total thruster flow equal the nominal demand
    (pump flow = thruster flow + orifice back flow).  Root found with the Illinois
    variant of regula falsi on log(diameter)   [DOWELL].  Returns (diameter or None, info)."""
    lower_evaluation = evaluate_orifice(model, orifice_identifier, lower_bound)
    if math.isnan(lower_evaluation['relative_error']):
        return None, {'message': 'No thrusters are ON in this scenario.', 'lower': lower_evaluation}
    if lower_evaluation['relative_error'] < 0:
        return None, {'message': 'Thruster demand cannot be met even with the orifice almost closed '
                                 '(pump flow too low or line losses too high).', 'lower': lower_evaluation}
    upper_evaluation = evaluate_orifice(model, orifice_identifier, upper_bound)
    if upper_evaluation['relative_error'] > 0:
        return None, {'message': f'Thruster flow still above nominal with a {upper_bound}" orifice.',
                      'upper': upper_evaluation}
    low, high = math.log(lower_bound), math.log(upper_bound)
    low_value, high_value = lower_evaluation['relative_error'], upper_evaluation['relative_error']
    guess = lower_evaluation['results']['initial_guess']
    last_side = 0
    evaluation = None
    for _ in range(80):
        # regula falsi point   [DOWELL]
        trial = (low * high_value - high * low_value) / (high_value - low_value) if high_value != low_value \
            else 0.5 * (low + high)
        if not (min(low, high) < trial < max(low, high)):
            trial = 0.5 * (low + high)
        evaluation = evaluate_orifice(model, orifice_identifier, math.exp(trial), guess)
        guess = evaluation['results']['initial_guess']
        trial_value = evaluation['relative_error']
        if abs(trial_value) < 1e-7 or abs(high - low) < 1e-10:
            break
        if trial_value > 0:
            low, low_value = trial, trial_value
            if last_side == -1:
                high_value *= 0.5          # Illinois modification   [DOWELL]
            last_side = -1
        else:
            high, high_value = trial, trial_value
            if last_side == 1:
                low_value *= 0.5           # Illinois modification   [DOWELL]
            last_side = 1
    return evaluation['diameter'], {'message': 'ok', 'evaluation': evaluation}


def size_orifice_for_all_scenarios(model, orifice_identifier, snapshots, tolerance_percent=2.0,
                                   grid_point_count=61, progress=None):
    """Find a single orifice diameter for several scenarios.

    Objective (lexicographic): fewest scenarios with limit violations, then the
    smallest worst-case |thruster flow / nominal - 1|.  A logarithmic grid scan is
    refined with golden-section search   [KIEFER]; the acceptable window edges are
    refined by bisection.  Also reports each scenario's own ideal diameter.
    """
    tolerance = tolerance_percent / 100.0
    scenario_models = []
    for snapshot in snapshots:
        scenario_model = apply_snapshot(copy.deepcopy(model), snapshot)
        density = float(scenario_model.get('fluid', DEFAULT_FLUID)['density'])
        if thruster_demand(scenario_model, density)[0] > 0:
            scenario_models.append((snapshot['name'], scenario_model))
    if not scenario_models:
        return {'message': 'No scenario has a thruster ON.', 'rows': [], 'best_diameter': None}
    guesses = [None] * len(scenario_models)

    def objective(diameter):
        evaluations = []
        for scenario_index, (_, scenario_model) in enumerate(scenario_models):
            evaluation = evaluate_orifice(scenario_model, orifice_identifier, diameter, guesses[scenario_index])
            guesses[scenario_index] = evaluation['results']['initial_guess']
            evaluations.append(evaluation)
        violating_count = sum(1 for evaluation in evaluations if evaluation['violations'])
        worst_error = max(abs(evaluation['relative_error']) for evaluation in evaluations)
        return violating_count * 1e3 + worst_error, evaluations

    # coarse logarithmic grid, then golden-section refinement in log(diameter)        [KIEFER]
    grid = numpy.exp(numpy.linspace(math.log(ORIFICE_DIAMETER_LOWER_BOUND),
                                    math.log(ORIFICE_DIAMETER_UPPER_BOUND), grid_point_count))
    grid_values = []
    for grid_index, diameter in enumerate(grid):
        grid_values.append(objective(diameter))
        if progress:
            progress(int(70 * (grid_index + 1) / grid_point_count))
    best_index = int(numpy.argmin([value[0] for value in grid_values]))
    low = math.log(grid[max(best_index - 1, 0)])
    high = math.log(grid[min(best_index + 1, grid_point_count - 1)])
    golden_ratio_fraction = (math.sqrt(5) - 1) / 2                     # [KIEFER]
    inner_low = high - golden_ratio_fraction * (high - low)
    inner_high = low + golden_ratio_fraction * (high - low)
    inner_low_value = objective(math.exp(inner_low))[0]
    inner_high_value = objective(math.exp(inner_high))[0]
    for refinement_index in range(40):
        if inner_low_value < inner_high_value:
            high, inner_high, inner_high_value = inner_high, inner_low, inner_low_value
            inner_low = high - golden_ratio_fraction * (high - low)
            inner_low_value = objective(math.exp(inner_low))[0]
        else:
            low, inner_low, inner_low_value = inner_low, inner_high, inner_high_value
            inner_high = low + golden_ratio_fraction * (high - low)
            inner_high_value = objective(math.exp(inner_high))[0]
        if progress:
            progress(70 + int(20 * (refinement_index + 1) / 40))
    best_diameter = math.exp(0.5 * (low + high))
    best_objective, best_evaluations = objective(best_diameter)

    def acceptable(diameter):
        _, evaluations = objective(diameter)
        return all(abs(evaluation['relative_error']) <= tolerance and not evaluation['violations']
                   for evaluation in evaluations)

    def bisect_edge(good_diameter, bad_diameter):
        for _ in range(20):
            middle = math.sqrt(good_diameter * bad_diameter)        # bisection in log(diameter) [DENNIS]
            if acceptable(middle):
                good_diameter = middle
            else:
                bad_diameter = middle
        return good_diameter

    window = None
    if acceptable(best_diameter):
        index = best_index
        while index > 0 and acceptable(grid[index - 1]):
            index -= 1
        lower_unacceptable = grid[index - 1] if index > 0 else None
        window_low = grid[index] if grid[index] < best_diameter else best_diameter
        index = best_index
        while index < grid_point_count - 1 and acceptable(grid[index + 1]):
            index += 1
        upper_unacceptable = grid[index + 1] if index < grid_point_count - 1 else None
        window_high = grid[index] if grid[index] > best_diameter else best_diameter
        if lower_unacceptable is not None:
            window_low = bisect_edge(window_low, lower_unacceptable)
        if upper_unacceptable is not None:
            window_high = bisect_edge(window_high, upper_unacceptable)
        window = (window_low, window_high)
    rows = []
    for (name, scenario_model), evaluation in zip(scenario_models, best_evaluations):
        ideal_diameter, information = ideal_orifice_for_model(scenario_model, orifice_identifier)
        rows.append({'name': name, 'ideal_diameter': ideal_diameter, 'ideal_message': information['message'],
                     'evaluation': evaluation})
    if progress:
        progress(100)
    return {'message': 'ok', 'best_diameter': best_diameter, 'best_objective': best_objective,
            'window': window, 'rows': rows, 'tolerance_percent': tolerance_percent}


def find_lowest_pump_speed(model, pump_identifier, snapshots, tolerance_percent=2.0):
    """For each scenario, solve with the pump ON in curve mode at every curve speed
    (ascending).  'lowest' is the lowest speed at which every ON thruster gets at
    least (1 - tolerance) of its nominal flow; 'lowest_without_warnings' also
    requires no pressure-limit warnings anywhere in the system."""
    tolerance = tolerance_percent / 100.0
    pump_block = find_block(model, pump_identifier)
    speeds = sorted(float(key) for key in (pump_block.get('pump_curves') or {}))
    scenario_results = []
    for snapshot in snapshots:
        scenario_model = (apply_snapshot(copy.deepcopy(model), snapshot) if snapshot
                          else copy.deepcopy(model))
        density = float(scenario_model.get('fluid', DEFAULT_FLUID)['density'])
        rows, lowest, lowest_without_warnings, guess = [], None, None, None
        demand, per_thruster = thruster_demand(scenario_model, density)
        for speed in speeds:
            trial_model = copy.deepcopy(scenario_model)
            pump = find_block(trial_model, pump_identifier)
            pump['state'] = 0
            pump['properties']['revolutions_per_minute'] = speed
            results = solve_model(trial_model, initial_guess=guess)
            guess = results['initial_guess']
            flow_ratios = [results['blocks'][identifier].get('mass_flow', 0.0) / nominal
                           for identifier, nominal in per_thruster if nominal > 0]
            thruster_flow = sum(results['blocks'][identifier].get('mass_flow', 0.0)
                                for identifier, _ in per_thruster)
            violations = limit_violations(results)
            meets_demand = bool(flow_ratios) and min(flow_ratios) >= 1 - tolerance and bool(results['converged'])
            rows.append({'speed': speed, 'thruster_flow': float(thruster_flow), 'demand': demand,
                         'minimum_flow_ratio': float(min(flow_ratios)) if flow_ratios else float('nan'),
                         'pump_flow': float(results['blocks'][pump_identifier].get('mass_flow', 0.0)),
                         'pump_pressure_rise': float(results['blocks'][pump_identifier]['calculated'].get(
                             'calculated_pressure_rise', float('nan'))),
                         'violations': violations, 'meets_demand': meets_demand})
            if meets_demand and lowest is None:
                lowest = speed
            if meets_demand and not violations and lowest_without_warnings is None:
                lowest_without_warnings = speed
        scenario_results.append({'name': snapshot['name'] if snapshot else 'Current', 'rows': rows,
                                 'lowest': lowest, 'lowest_without_warnings': lowest_without_warnings,
                                 'demand': demand})
    return scenario_results


# =============================================================================
# Orthogonal polyline geometry and routing (pure Python, points are (x, y))
# =============================================================================
GRID = 10
GEOMETRY_TOLERANCE = 1e-6


def snap_to_grid(value, grid=GRID):
    return round(value / grid) * grid


def _raw_orientation(first_point, second_point):
    """'H' horizontal, 'V' vertical, 'D' diagonal, None for zero length."""
    if (abs(first_point[0] - second_point[0]) < GEOMETRY_TOLERANCE
            and abs(first_point[1] - second_point[1]) < GEOMETRY_TOLERANCE):
        return None
    if abs(first_point[1] - second_point[1]) < GEOMETRY_TOLERANCE:
        return 'H'
    if abs(first_point[0] - second_point[0]) < GEOMETRY_TOLERANCE:
        return 'V'
    return 'D'


def _flip_orientation(orientation):
    return 'V' if orientation == 'H' else 'H'


def segment_orientation(points, segment_index):
    """Orientation of a segment ('H' or 'V'); zero-length segments take the
    orientation implied by the alternating H/V pattern of their neighbours."""
    orientation = _raw_orientation(points[segment_index], points[segment_index + 1])
    if orientation is not None:
        return 'H' if orientation == 'D' else orientation
    for other_index in range(segment_index - 1, -1, -1):
        other = _raw_orientation(points[other_index], points[other_index + 1])
        if other in ('H', 'V'):
            return other if (segment_index - other_index) % 2 == 0 else _flip_orientation(other)
    for other_index in range(segment_index + 1, len(points) - 1):
        other = _raw_orientation(points[other_index], points[other_index + 1])
        if other in ('H', 'V'):
            return other if (other_index - segment_index) % 2 == 0 else _flip_orientation(other)
    return 'H'


def orthogonalize(points):
    """Insert a corner wherever two consecutive points are diagonal."""
    result = [tuple(points[0])]
    for point in points[1:]:
        previous_point = result[-1]
        if _raw_orientation(previous_point, point) == 'D':
            previous_orientation = (_raw_orientation(result[-2], previous_point) if len(result) >= 2 else None)
            corner = ((previous_point[0], point[1]) if previous_orientation == 'H'
                      else (point[0], previous_point[1]))
            result.append(corner)
        result.append(tuple(point))
    return result


def tidy(points):
    """Merge collinear neighbouring segments; keeps zero-length split markers."""
    points = [tuple(point) for point in points]
    changed = True
    while changed and len(points) > 2:
        changed = False
        for index in range(1, len(points) - 1):
            before = _raw_orientation(points[index - 1], points[index])
            after = _raw_orientation(points[index], points[index + 1])
            if before in ('H', 'V') and before == after:
                del points[index]
                changed = True
                break
    return points


def simplify(points):
    points = orthogonalize(points)
    result = [points[0]]
    for point in points[1:]:
        if _raw_orientation(result[-1], point) is not None:
            result.append(point)
    if len(result) == 1:
        result.append(tuple(points[-1]))
    return tidy(result)


def polyline_length(points):
    return sum(abs(points[index + 1][0] - points[index][0]) + abs(points[index + 1][1] - points[index][1])
               for index in range(len(points) - 1))


def point_at_fraction(points, path_fraction):
    """Point at a fraction of the path length, and the index of its segment."""
    total_length = polyline_length(points)
    if total_length <= 0:
        return tuple(points[0]), 0
    target, travelled = min(max(path_fraction, 0.0), 1.0) * total_length, 0.0
    for index in range(len(points) - 1):
        start, end = points[index], points[index + 1]
        segment_length = abs(end[0] - start[0]) + abs(end[1] - start[1])
        if segment_length > 0 and travelled + segment_length >= target - 1e-9:
            fraction = (target - travelled) / segment_length
            return (start[0] + (end[0] - start[0]) * fraction, start[1] + (end[1] - start[1]) * fraction), index
        travelled += segment_length
    return tuple(points[-1]), len(points) - 2


def project_on_polyline(points, point):
    """Closest point on a polyline -> (path fraction, point, segment index, distance)."""
    total_length = polyline_length(points)
    best = (0.0, tuple(points[0]), 0, float('inf'))
    travelled = 0.0
    for index in range(len(points) - 1):
        start, end = points[index], points[index + 1]
        delta_x, delta_y = end[0] - start[0], end[1] - start[1]
        squared_length = delta_x * delta_x + delta_y * delta_y
        segment_length = abs(delta_x) + abs(delta_y)
        along = (0.0 if squared_length == 0 else
                 max(0.0, min(1.0, ((point[0] - start[0]) * delta_x + (point[1] - start[1]) * delta_y)
                              / squared_length)))
        closest = (start[0] + along * delta_x, start[1] + along * delta_y)
        distance = math.hypot(point[0] - closest[0], point[1] - closest[1])
        if distance < best[3]:
            best = ((travelled + along * segment_length) / total_length if total_length > 0 else 0.0,
                    closest, index, distance)
        travelled += segment_length
    return best


def sub_polyline(points, start_fraction, end_fraction):
    start_point, start_index = point_at_fraction(points, start_fraction)
    end_point, end_index = point_at_fraction(points, end_fraction)
    result = ([start_point] + [tuple(points[index]) for index in range(start_index + 1, end_index + 1)]
              + [end_point])
    cleaned = [result[0]]
    for point in result[1:]:
        if _raw_orientation(cleaned[-1], point) is not None:
            cleaned.append(point)
    return cleaned if len(cleaned) > 1 else [start_point, end_point]


def move_segment(points, segment_index, offset, start_fixed=True, end_fixed=True):
    """Drag a segment perpendicular to itself by `offset`.  Fixed end points get a
    new stub segment so they stay put.  Returns (new points, new segment index)."""
    points = [tuple(point) for point in points]
    orientation = segment_orientation(points, segment_index)
    if segment_index == 0 and start_fixed:
        points.insert(0, points[0])
        segment_index += 1
    if segment_index == len(points) - 2 and end_fixed:
        points.append(points[-1])
    for index in (segment_index, segment_index + 1):
        x, y = points[index]
        points[index] = (x, y + offset) if orientation == 'H' else (x + offset, y)
    return points, segment_index


def adjust_endpoint(points, which_end, new_position):
    """Move one end point, keeping the polyline orthogonal."""
    points = [tuple(point) for point in points]
    reversed_order = which_end == 'start'
    if reversed_order:
        points.reverse()
    new_position = tuple(new_position)
    if len(points) == 2:
        fixed_point = points[0]
        orientation = segment_orientation(points, 0)
        corner = ((new_position[0], fixed_point[1]) if orientation == 'H'
                  else (fixed_point[0], new_position[1]))
        points = [fixed_point, corner, new_position]
    else:
        orientation = segment_orientation(points, len(points) - 2)
        previous_x, previous_y = points[-2]
        points[-2] = (previous_x, new_position[1]) if orientation == 'H' else (new_position[0], previous_y)
        points[-1] = new_position
    if reversed_order:
        points.reverse()
    return orthogonalize(points)


def orthogonal_connection(last_point, target, previous_orientation=None):
    """Points to append to go orthogonally from `last_point` to `target`."""
    if _raw_orientation(last_point, target) in (None, 'H', 'V'):
        return [tuple(target)]
    corner = ((last_point[0], target[1]) if previous_orientation == 'H'
              else (target[0], last_point[1]))
    return [corner, tuple(target)]


def split_segment(points, segment_index, position):
    """'Add segment': split a segment at `position` (a zero-length joint that can
    then be dragged out)."""
    points = [tuple(point) for point in points]
    start, end = points[segment_index], points[segment_index + 1]
    if segment_orientation(points, segment_index) == 'H':
        x = min(max(snap_to_grid(position[0]), min(start[0], end[0])), max(start[0], end[0]))
        split_point = (x, start[1])
    else:
        y = min(max(snap_to_grid(position[1]), min(start[1], end[1])), max(start[1], end[1]))
        split_point = (start[0], y)
    return points[:segment_index + 1] + [split_point, split_point] + points[segment_index + 1:]


def delete_segment(points, segment_index, start_fixed, end_fixed):
    """Remove a segment by sliding a neighbour; returns new points or None."""
    points = [tuple(point) for point in points]
    segment_count = len(points) - 1
    if segment_count <= 1:
        return None
    if segment_index == 0 and not start_fixed:
        return simplify(points[1:])
    if segment_index == segment_count - 1 and not end_fixed:
        return simplify(points[:-1])
    orientation = segment_orientation(points, segment_index)
    start, end = points[segment_index], points[segment_index + 1]
    length_along = (end[0] - start[0]) if orientation == 'H' else (end[1] - start[1])

    def shift(point_list, indices, distance):
        for index in indices:
            x, y = point_list[index]
            point_list[index] = (x + distance, y) if orientation == 'H' else (x, y + distance)
    if segment_index + 1 <= segment_count - 1 and not (segment_index + 1 == segment_count - 1 and end_fixed):
        new_points = list(points)
        shift(new_points, (segment_index + 1, segment_index + 2), -length_along)
        return simplify(new_points)
    if segment_index - 1 >= 0 and not (segment_index - 1 == 0 and start_fixed):
        new_points = list(points)
        shift(new_points, (segment_index - 1, segment_index), length_along)
        return simplify(new_points)
    return None


def path_crosses_obstacles(points, rectangles, port_positions=(), shrink=2.0, port_radius=3.0):
    """True if any segment passes over a block body or over a foreign port."""
    for index in range(len(points) - 1):
        start, end = points[index], points[index + 1]
        left, right = sorted((start[0], end[0]))
        top, bottom = sorted((start[1], end[1]))
        for rectangle_left, rectangle_top, rectangle_right, rectangle_bottom in rectangles:
            if (right > rectangle_left + shrink and left < rectangle_right - shrink
                    and bottom > rectangle_top + shrink and top < rectangle_bottom - shrink):
                return True
        for port_position in port_positions:
            if (left - port_radius <= port_position[0] <= right + port_radius
                    and top - port_radius <= port_position[1] <= bottom + port_radius):
                return True
    return False


GRID_DIRECTIONS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def route_around_blocks(start, start_direction, end, end_direction, rectangles, blocked_positions=(),
                        step=GRID, stub_cells=2, clearance=10, bend_penalty=4.0, maximum_cells=400000):
    """Orthogonal A* route [HART] on the grid that never crosses a block body
    (inflated by `clearance`) or a foreign port.  Cost = path length + a penalty
    per bend; heuristic = Manhattan distance (admissible).  Start and end
    directions are outward port directions or None.  Returns points or None."""
    def cell_of(point):
        return (int(round(point[0] / step)), int(round(point[1] / step)))
    start_cell, end_cell = cell_of(start), cell_of(end)
    route_start = ((start_cell[0] + start_direction[0] * stub_cells, start_cell[1] + start_direction[1] * stub_cells)
                   if start_direction else start_cell)
    route_end = ((end_cell[0] + end_direction[0] * stub_cells, end_cell[1] + end_direction[1] * stub_cells)
                 if end_direction else end_cell)
    blocked = set()
    x_values = [route_start[0], route_end[0], start_cell[0], end_cell[0]]
    y_values = [route_start[1], route_end[1], start_cell[1], end_cell[1]]
    for left, top, right, bottom in rectangles:
        first_column = math.floor((left - clearance) / step + 1e-9) + 1
        last_column = math.ceil((right + clearance) / step - 1e-9) - 1
        first_row = math.floor((top - clearance) / step + 1e-9) + 1
        last_row = math.ceil((bottom + clearance) / step - 1e-9) - 1
        x_values += [first_column, last_column]
        y_values += [first_row, last_row]
        for column in range(first_column, last_column + 1):
            for row in range(first_row, last_row + 1):
                blocked.add((column, row))
    for position in blocked_positions:
        blocked.add(cell_of(position))
    blocked.discard(route_start)
    blocked.discard(route_end)
    minimum_column, maximum_column = min(x_values) - 8, max(x_values) + 8
    minimum_row, maximum_row = min(y_values) - 8, max(y_values) + 8
    if (maximum_column - minimum_column) * (maximum_row - minimum_row) > maximum_cells:
        return None
    start_heading = GRID_DIRECTIONS.index(tuple(start_direction)) if start_direction else -1
    goal_heading = GRID_DIRECTIONS.index((-end_direction[0], -end_direction[1])) if end_direction else -1

    def heuristic(cell):
        return abs(cell[0] - route_end[0]) + abs(cell[1] - route_end[1])
    initial_state = (route_start[0], route_start[1], start_heading)
    cost_so_far = {initial_state: 0.0}
    parent = {initial_state: None}
    frontier = [(heuristic(route_start), 0.0, initial_state)]
    found = None
    while frontier:
        _, current_cost, state = heapq.heappop(frontier)
        if current_cost > cost_so_far.get(state, float('inf')):
            continue
        column, row, heading = state
        if (column, row) == route_end:
            found = state
            break
        for new_heading, (step_x, step_y) in enumerate(GRID_DIRECTIONS):
            if heading != -1 and new_heading == (heading + 2) % 4:
                continue   # no reversing
            next_column, next_row = column + step_x, row + step_y
            if (next_column < minimum_column or next_column > maximum_column or next_row < minimum_row
                    or next_row > maximum_row or (next_column, next_row) in blocked):
                continue
            new_cost = current_cost + 1.0 + (bend_penalty if (heading != -1 and new_heading != heading) else 0.0)
            if (next_column, next_row) == route_end and goal_heading != -1 and new_heading != goal_heading:
                new_cost += bend_penalty
            next_state = (next_column, next_row, new_heading)
            if new_cost < cost_so_far.get(next_state, float('inf')):
                cost_so_far[next_state] = new_cost
                parent[next_state] = state
                heapq.heappush(frontier, (new_cost + heuristic((next_column, next_row)), new_cost, next_state))
    if found is None:
        return None
    cells = []
    state = found
    while state is not None:
        cells.append((state[0] * step, state[1] * step))
        state = parent[state]
    cells.reverse()
    return simplify([tuple(start)] + cells + [tuple(end)])


def quick_orthogonal_route(start, start_direction, end, end_direction, stub_length=20):
    """Quick (non-avoiding) orthogonal route used for live previews."""
    start_stub = ((start[0] + start_direction[0] * stub_length, start[1] + start_direction[1] * stub_length)
                  if start_direction else tuple(start))
    end_stub = ((end[0] + end_direction[0] * stub_length, end[1] + end_direction[1] * stub_length)
                if end_direction else tuple(end))
    return simplify([tuple(start), start_stub, (end_stub[0], start_stub[1]), end_stub, tuple(end)])


def _overlap_area(first, second):
    width = min(first[2], second[2]) - max(first[0], second[0])
    height = min(first[3], second[3]) - max(first[1], second[1])
    return width * height if width > 0 and height > 0 else 0.0


def place_label(segment_start, segment_end, width, height, obstacles, gap=9):
    """Centre for a label beside a segment that avoids `obstacles`
    (left, top, right, bottom).  Tries above/below (or right/left) at several
    positions along the segment and returns the first clear spot, else the least overlap."""
    horizontal = abs(segment_start[1] - segment_end[1]) < 1e-6
    middle_x, middle_y = (segment_start[0] + segment_end[0]) / 2, (segment_start[1] + segment_end[1]) / 2
    length = abs(segment_end[0] - segment_start[0]) + abs(segment_end[1] - segment_start[1])
    candidates = []
    for shift in (0.0, -0.25 * length, 0.25 * length, -0.42 * length, 0.42 * length):
        if horizontal:
            candidates += [(middle_x + shift, middle_y - gap - height / 2),
                           (middle_x + shift, middle_y + gap + height / 2)]
        else:
            candidates += [(middle_x + gap + width / 2, middle_y + shift),
                           (middle_x - gap - width / 2, middle_y + shift)]
    best = None
    for candidate_index, (center_x, center_y) in enumerate(candidates):
        rectangle = (center_x - width / 2, center_y - height / 2, center_x + width / 2, center_y + height / 2)
        overlap = sum(_overlap_area(rectangle, obstacle) for obstacle in obstacles)
        if overlap == 0:
            return center_x, center_y
        if best is None or overlap + candidate_index < best[0]:
            best = (overlap + candidate_index, (center_x, center_y))
    return best[1]


# ---- block geometry (shared by the GUI and the headless example builder) ----
def block_size(block_type, orientation=0):
    width, height = BLOCK_TYPES[block_type]['size']
    return (height, width) if orientation % 180 == 90 else (width, height)


def rotate_vector(x, y, orientation):
    """Rotate by 90° steps, clockwise on screen (y axis pointing down)."""
    angle = orientation % 360
    if angle == 90:
        return -y, x
    if angle == 180:
        return -x, -y
    if angle == 270:
        return y, -x
    return x, y


def port_local_position(block_type, port_name, orientation=0):
    x, y, direction_x, direction_y = BLOCK_TYPES[block_type]['ports'][port_name]
    position_x, position_y = rotate_vector(x, y, orientation)
    rotated_direction_x, rotated_direction_y = rotate_vector(direction_x, direction_y, orientation)
    return (position_x, position_y), (int(round(rotated_direction_x)), int(round(rotated_direction_y)))


def block_rectangle(block):
    width, height = block_size(block['type'], block.get('orientation', 0))
    return (block['x'] - width / 2, block['y'] - height / 2, block['x'] + width / 2, block['y'] + height / 2)


def port_scene_position(block, port_name):
    (position_x, position_y), direction = port_local_position(block['type'], port_name,
                                                              block.get('orientation', 0))
    return (block['x'] + position_x, block['y'] + position_y), direction


def route_model_lines(model):
    """Give every line without points an A* route (used for generated models),
    and resolve taps given as an 'at' point into a path fraction."""
    blocks = blocks_by_identifier(model)
    lines = {line['identifier']: line for line in model.get('lines', [])}
    rectangles = [block_rectangle(block) for block in blocks.values()]
    all_ports = [port_scene_position(block, port_name)[0] for block in blocks.values()
                 for port_name in BLOCK_TYPES[block['type']]['ports']]

    def endpoint(attachment):
        if attachment and attachment.get('kind') == 'port' and attachment.get('block') in blocks:
            return port_scene_position(blocks[attachment['block']], attachment['port'])
        if attachment and attachment.get('kind') == 'line' and 'at' in attachment:
            return tuple(attachment['at']), None
        if (attachment and attachment.get('kind') == 'line' and attachment.get('line') in lines
                and lines[attachment['line']].get('points')):
            return point_at_fraction(lines[attachment['line']]['points'],
                                     float(attachment.get('path_fraction', 0)))[0], None
        return None, None
    order = sorted(lines.values(), key=lambda line: sum(1 for attachment in (line.get('start'), line.get('end'))
                                                        if attachment and attachment.get('kind') == 'line'))
    for line in order:
        if line.get('points'):
            continue
        start, start_direction = endpoint(line.get('start'))
        end, end_direction = endpoint(line.get('end'))
        if start is None or end is None:
            continue
        foreign_ports = [port for port in all_ports if port != tuple(start) and port != tuple(end)]
        points = (route_around_blocks(start, start_direction, end, end_direction, rectangles, foreign_ports)
                  or quick_orthogonal_route(start, start_direction, end, end_direction))
        line['points'] = [list(point) for point in points]
    for line in lines.values():
        for end_name in ('start', 'end'):
            attachment = line.get(end_name)
            if attachment and attachment.get('kind') == 'line' and 'at' in attachment:
                target = lines.get(attachment['line'])
                if target and target.get('points'):
                    attachment['path_fraction'] = project_on_polyline(
                        [tuple(point) for point in target['points']], tuple(attachment['at']))[0]
                attachment.pop('at', None)
    return model


# =============================================================================
# Example system
# =============================================================================
def incremented_name(name, taken_names):
    """Name for a copy of a block: the last number in the name plus one ('LV-1' -> 'LV-2',
    'TK-09' -> 'TK-10', 'Valve 2 left' -> 'Valve 3 left'); a name without a number gets ' 1'
    ('Main valve' -> 'Main valve 1').  Keeps counting past names already in use."""
    match = None
    for match in re.finditer(r'\d+', name):
        pass                                             # keep the last run of digits
    if match is not None:
        before, after = name[:match.start()], name[match.end():]
        number, width = int(match.group()) + 1, len(match.group())
    else:
        before = name if (not name or name[-1] in ' -_') else name + ' '
        after, number, width = '', 1, 1
    while True:
        candidate = f"{before}{number:0{width}d}{after}"
        if candidate not in taken_names:
            return candidate
        number += 1


def new_identifier():
    return uuid.uuid4().hex[:8]


def build_example_model():
    """Pump-fed example: tank, filter, latch valve, pump with a bypass orifice
    back to the pump inlet manifold, flowmeter and thruster."""
    def make_block(identifier, block_type, name, x, y, orientation=0, **property_overrides):
        properties = default_block_properties(block_type)
        properties.update(property_overrides)
        update_linked_flow_coefficients(block_type, properties)
        return {'identifier': identifier, 'type': block_type, 'name': name, 'x': x, 'y': y,
                'orientation': orientation, 'state': 0, 'properties': properties, 'pump_curves': {}}
    quarter_inch_tube = '1/4" x 0.035"'
    eighth_inch_tube = '1/8" x 0.028"'
    blocks = [
        make_block('tk1', 'tank', 'TK-1', 120, -20, volume=20.0),
        make_block('fl1', 'filter', 'FL-1', 340, 300),
        make_block('lv1', 'latch_valve', 'LV-1', 600, 300),
        make_block('t1', 'tee', 'T-1', 860, 300),
        make_block('pt1', 'transducer', 'PT-1', 740, 150),
        make_block('p1', 'pump', 'P-1', 1100, 300),
        make_block('pt2', 'transducer', 'PT-2', 1260, 150),
        make_block('t2', 'tee', 'T-2', 1400, 300),
        make_block('fm1', 'flowmeter', 'FM-1', 1640, 300, outlet_tubing_size=eighth_inch_tube),
        make_block('th1', 'thruster', 'TH-1', 1900, 300, inlet_tubing_size=eighth_inch_tube),
        make_block('or1', 'orifice', 'OR-1', 1130, 520, orientation=180, orifice_diameter=0.00914),
    ]
    # synthetic gear-pump curves: Q = 0.12 mL/rev * rpm - 0.12 mL/min/psid * pressure rise
    pump_curves = {}
    for speed in (2000, 3000, 4000, 5000, 6000):
        pump_curves[str(speed)] = [[pressure_rise, round(0.12 * speed - 0.12 * pressure_rise, 3)]
                                   for pressure_rise in range(0, 601, 100)]
    blocks[5]['pump_curves'] = pump_curves

    def make_line(identifier, name, start, end, length=6.0, bend_count=0, tube=quarter_inch_tube):
        outer_diameter, wall_thickness = tube_outer_diameter_and_wall(tube)
        properties = default_line_properties()
        properties.update(outer_diameter=outer_diameter, inner_diameter=outer_diameter - 2 * wall_thickness,
                          length=length, bend_count=bend_count)
        return {'identifier': identifier, 'name': name, 'points': None, 'start': start, 'end': end,
                'properties': properties, 'auto_routed': True}

    def port(block_identifier, port_name):
        return {'kind': 'port', 'block': block_identifier, 'port': port_name}
    lines = [
        make_line('l1', 'L1', port('tk1', 'outlet'), port('fl1', 'inlet'), 24.0, 2),
        make_line('l2', 'L2', port('fl1', 'outlet'), port('lv1', 'inlet'), 6.0),
        make_line('l3', 'L3', port('lv1', 'outlet'), port('t1', 'inlet'), 6.0, 1),
        make_line('l4', 'L4', port('t1', 'outlet'), port('p1', 'inlet'), 6.0),
        make_line('l5', 'L5', port('p1', 'outlet'), port('t2', 'inlet'), 6.0),
        make_line('l6', 'L6', port('t2', 'outlet'), port('fm1', 'inlet'), 12.0, 1),
        make_line('l7', 'L7', port('fm1', 'outlet'), port('th1', 'inlet'), 18.0, 3, eighth_inch_tube),
        make_line('l8', 'L8', port('t2', 'branch'), port('or1', 'inlet'), 8.0, 2),
        make_line('l9', 'L9', port('or1', 'outlet'), port('t1', 'branch'), 8.0, 2),
        make_line('l10', 'L10', port('pt1', 'inlet'), {'kind': 'line', 'line': 'l3', 'at': [740, 300]}, 3.0, 1),
        make_line('l11', 'L11', port('pt2', 'inlet'), {'kind': 'line', 'line': 'l5', 'at': [1260, 300]}, 3.0, 1),
    ]
    model = {'format': 'pylon_schematic', 'version': FILE_FORMAT_VERSION, 'fluid': dict(DEFAULT_FLUID),
             'settings': dict(DEFAULT_SETTINGS), 'blocks': blocks, 'lines': lines,
             'scenarios': [], 'active_scenario': None}
    route_model_lines(model)
    nominal = make_snapshot(model, 'Nominal')
    maximum_pressure_scenario = copy.deepcopy(nominal)
    maximum_pressure_scenario['name'] = 'Maximum tank pressure, clogged filter'
    maximum_pressure_scenario['blocks']['tk1']['state'] = 1
    maximum_pressure_scenario['blocks']['fl1']['state'] = 1
    minimum_pressure_scenario = copy.deepcopy(nominal)
    minimum_pressure_scenario['name'] = 'Minimum tank pressure'
    minimum_pressure_scenario['blocks']['tk1']['state'] = 2
    model['scenarios'] = [nominal, maximum_pressure_scenario, minimum_pressure_scenario]
    model['active_scenario'] = 'Nominal'
    return model


# =============================================================================
# Self test (python pylon.py --selftest ; no GUI needed)
# =============================================================================
def run_self_test():
    all_passed = True

    def check(name, actual, expected, relative_tolerance):
        nonlocal all_passed
        passed = abs(actual - expected) <= relative_tolerance * abs(expected)
        all_passed &= passed
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}: got {actual:.6g}, expected {expected:.6g}")
    print("Physics checks")
    density, viscosity = 998.2, 1.002e-3          # water at 20 C   [WHITE]
    # 1. thruster equivalent orifice against the closed form  mdot = Cd A sqrt(2 rho dP)   [SUTTON]
    model = {'fluid': {'name': 'Water', 'density': density, 'viscosity_millipascal_seconds': 1.002},
             'settings': {},
             'blocks': [{'identifier': 'a', 'type': 'tank', 'state': 0,
                         'properties': default_block_properties('tank')},
                        {'identifier': 'b', 'type': 'thruster', 'state': 0,
                         'properties': default_block_properties('thruster')}],
             'lines': [{'identifier': 'L', 'start': {'kind': 'port', 'block': 'a', 'port': 'outlet'},
                        'end': {'kind': 'port', 'block': 'b', 'port': 'inlet'},
                        'properties': dict(default_line_properties(), length=1e-6)}]}
    model['blocks'][0]['properties'].update(nominal_pressure=200.0, outlet_tubing_size=DEFAULT_TUBE)
    model['blocks'][1]['properties'].update(orifice_diameter=0.02, inlet_tubing_size=DEFAULT_TUBE)
    results = solve_model(model)
    orifice_area = math.pi * (0.02 * METERS_PER_INCH) ** 2 / 4
    # the tank entrance loss on the tube area is negligible next to the orifice
    check('Equivalent-orifice thruster flow', results['blocks']['b']['mass_flow'],
          0.6 * orifice_area * math.sqrt(2 * density * 200 * PASCALS_PER_PSI), 2e-3)
    # 2. laminar tube against Hagen-Poiseuille  dP = 128 mu L Q / (pi D^4)   [WHITE]
    diameter, length, mass_flow = 0.18 * METERS_PER_INCH, 100 * METERS_PER_INCH, 1e-4
    pressure_drop_function = make_tube_pressure_drop_function(diameter, length, 0, 3, 0.0, density, viscosity,
                                                              DEFAULT_SETTINGS)
    volumetric_flow = mass_flow / density
    check('Laminar tube pressure drop (Hagen-Poiseuille)', pressure_drop_function(mass_flow),
          128 * viscosity * length * volumetric_flow / (math.pi * diameter ** 4), 1e-6)
    # 3. turbulent tube against Colebrook  1/sqrt(f) = -2 log10(eps/3.7D + 2.51/(Re sqrt f))   [COLEBROOK]
    mass_flow = 0.2
    area = math.pi * diameter * diameter / 4
    velocity = mass_flow / (density * area)
    reynolds_number = density * velocity * diameter / viscosity
    friction_factor = 0.02
    for _ in range(50):
        friction_factor = (-2 * math.log10(1.5e-6 / diameter / 3.7
                                           + 2.51 / (reynolds_number * math.sqrt(friction_factor)))) ** -2
    check(f'Turbulent tube pressure drop vs Colebrook (Re={reynolds_number:.0f})',
          make_tube_pressure_drop_function(diameter, length, 0, 3, 1.5e-6, density, viscosity,
                                           DEFAULT_SETTINGS)(mass_flow),
          friction_factor * length / diameter * density * velocity * velocity / 2, 0.03)   # Darcy [CRANE]
    # 4. Cv relation   [ISA]
    check('Cv = 29.84 Cd d^2', flow_coefficient_from_orifice_diameter(0.1, 0.6), 29.84 * 0.6 * 0.01, 2e-3)
    print("Network checks (example system)")
    example = build_example_model()
    results = solve_model(example)
    print(f"  solver: {results['message']}")
    all_passed &= results['converged']
    block_results = results['blocks']
    pump_flow = block_results['p1']['mass_flow']
    thruster_flow = block_results['th1']['mass_flow']
    orifice_flow = block_results['or1']['mass_flow']
    check('Pump flow = thruster flow + orifice back flow', thruster_flow + orifice_flow, pump_flow, 1e-8)
    tee_calculated = block_results['t2']['calculated']
    check('Tee T-2 mass balance', tee_calculated['calculated_outlet_flow'] + tee_calculated['calculated_branch_flow'],
          tee_calculated['calculated_inlet_flow'], 1e-8)
    print(f"  PT-1 {block_results['pt1']['reading']:.2f} psia, PT-2 {block_results['pt2']['reading']:.2f} psia, "
          f"thruster {thruster_flow * 1e3:.3f} g/s, bypass {orifice_flow * 1e3:.3f} g/s")
    ideal_diameter, information = ideal_orifice_for_model(example, 'or1')
    all_passed &= ideal_diameter is not None
    if ideal_diameter:
        evaluation = information['evaluation']
        check(f'Ideal bypass orifice ({ideal_diameter:.5f} in) gives nominal thruster flow',
              evaluation['thruster_flow'], evaluation['demand'], 1e-5)
    print("Routing checks")
    route = route_around_blocks((0, 0), (1, 0), (300, 0), (-1, 0), [(100, -40, 200, 40)])
    crosses = path_crosses_obstacles(route, [(100, -40, 200, 40)])
    print(f"  route around block: {route}")
    all_passed &= (route is not None and not crosses)
    print("File upgrade check")
    old_model = {'version': 2, 'fluid': {'name': 'Water', 'rho': 998.2, 'mu_mPas': 1.002},
                 'settings': {'orifice_cd': 0.6},
                 'blocks': [{'id': 'x', 'type': 'pump', 'x': 0, 'y': 0, 'orient': 0, 'state': 0,
                             'props': {'flow': 5.0, 'eseod_stop': 0.02}, 'curves': {}}],
                 'lines': [{'id': 'y', 'points': [[0, 0], [10, 0]], 'auto': True,
                            'start': {'kind': 'port', 'block': 'x', 'port': 'out'}, 'end': None,
                            'props': {'od': 0.25, 'id': 0.18, 'bends': 1}}],
                 'scenarios': [{'name': 's', 'blocks': {'x': {'state': 0, 'props': {'flow': 4.0}}},
                                'lines': {'y': {'props': {'bends': 2}}}}]}
    upgraded = upgrade_model(old_model)
    upgrade_ok = (upgraded['blocks'][0]['properties'].get('stopped_orifice_diameter') == 0.02
                  and 'flow_rate' not in upgraded['blocks'][0]['properties']
                  and upgraded['lines'][0]['properties'].get('bend_count') == 1
                  and upgraded['lines'][0]['start']['port'] == 'outlet'
                  and 'flow_rate' not in upgraded['scenarios'][0]['blocks']['x']['properties']
                  and upgraded['fluid']['density'] == 998.2
                  and upgraded['settings'].get('orifice_discharge_coefficient') == 0.6)
    print(f"  [{'PASS' if upgrade_ok else 'FAIL'}] version 2 file upgraded to version {upgraded['version']}")
    all_passed &= upgrade_ok

    print("Pump curve between stored speeds")
    # gear-pump curves Q = 0.12 * N - 0.12 * dP (mL/min) stored at 3000 and 4000 rpm only;
    # at 3500 rpm and 200 psid the flow must be 0.12*3500 - 0.12*200 = 396 mL/min
    speed_curves = {str(speed): [[rise, 0.12 * speed - 0.12 * rise] for rise in range(0, 601, 100)]
                    for speed in (3000, 4000)}
    curve_function, _ = make_pump_curve_function_for_speed(speed_curves, 3500.0, 1000.0)
    interpolated_flow = curve_function(200 * PASCALS_PER_PSI)[0] / 1000.0 / CUBIC_METERS_PER_SECOND_PER_MILLILITER_PER_MINUTE
    check('Pump flow at 3500 rpm, 200 psid [mL/min]', interpolated_flow, 396.0, 1e-9)
    _, outside_note = make_pump_curve_function_for_speed(speed_curves, 5000.0, 1000.0)
    outside_ok = outside_note is not None
    print(f"  [{'PASS' if outside_ok else 'FAIL'}] speed outside the stored curves is reported")
    all_passed &= outside_ok

    print("Cross fitting")
    # tank -> cross inlet; identical thrusters on the outlet and both branches: the three
    # branches through the branch arms must carry equal flow and the cross must balance
    cross_model = {'fluid': dict(DEFAULT_FLUID), 'settings': {},
                   'blocks': [{'identifier': 'tank', 'type': 'tank', 'state': 0,
                               'properties': dict(default_block_properties('tank'), nominal_pressure=300.0)},
                              {'identifier': 'cross', 'type': 'cross', 'state': 0,
                               'properties': default_block_properties('cross')}]
                   + [{'identifier': name, 'type': 'thruster', 'state': 0,
                       'properties': dict(default_block_properties('thruster'), inlet_tubing_size=DEFAULT_TUBE)}
                      for name in ('east', 'north', 'south')],
                   'lines': []}
    for line_name, start, end in (('feed', ('tank', 'outlet'), ('cross', 'inlet')),
                                  ('east_line', ('cross', 'outlet'), ('east', 'inlet')),
                                  ('north_line', ('cross', 'upper_branch'), ('north', 'inlet')),
                                  ('south_line', ('cross', 'lower_branch'), ('south', 'inlet'))):
        cross_model['lines'].append({'identifier': line_name, 'points': [[0, 0], [10, 0]],
                                     'start': {'kind': 'port', 'block': start[0], 'port': start[1]},
                                     'end': {'kind': 'port', 'block': end[0], 'port': end[1]},
                                     'properties': dict(default_line_properties(), length=6.0)})
    cross_results = solve_model(cross_model)
    cross_calculated = cross_results['blocks']['cross']['calculated']
    check('Cross mass balance (in = out + up + down)',
          cross_calculated['calculated_outlet_flow'] + cross_calculated['calculated_upper_branch_flow']
          + cross_calculated['calculated_lower_branch_flow'], cross_calculated['calculated_inlet_flow'], 1e-8)
    check('Cross branches symmetric (up = down)', cross_calculated['calculated_upper_branch_flow'],
          cross_calculated['calculated_lower_branch_flow'], 1e-6)
    print("COPV, fluid coupling, fill & drain valve")

    def small_model(block_list, connections):
        built = {'fluid': dict(DEFAULT_FLUID), 'settings': {}, 'lines': [],
                 'blocks': [{'identifier': name, 'type': block_type, 'state': state,
                             'properties': dict(default_block_properties(block_type), **overrides)}
                            for name, block_type, state, overrides in block_list]}
        for index, (start, end) in enumerate(connections):
            built['lines'].append({'identifier': f"line{index}", 'points': [[0, 0], [10, 0]],
                                   'start': {'kind': 'port', 'block': start[0], 'port': start[1]},
                                   'end': {'kind': 'port', 'block': end[0], 'port': end[1]},
                                   'properties': dict(default_line_properties(), length=6.0)})
        return built
    feed = [(('tank', 'outlet'), ('engine', 'inlet'))]
    pressurant_blocks = [('tank', 'tank', 0, {'nominal_pressure': 50.0}),
                         ('engine', 'thruster', 0, {'inlet_tubing_size': DEFAULT_TUBE})]
    backfilled = solve_model(small_model(
        pressurant_blocks + [('bottle', 'copv', 0, {'ullage_pressure_setpoint': 300.0})],
        feed + [(('bottle', 'outlet'), ('tank', 'inlet'))]))
    check('COPV backfill sets tank ullage pressure [psia]',
          backfilled['blocks']['tank']['calculated']['calculated_pressure'], 300.0, 1e-12)
    isolated = solve_model(small_model(
        pressurant_blocks + [('bottle', 'copv', 1, {'ullage_pressure_setpoint': 300.0})],
        feed + [(('bottle', 'outlet'), ('tank', 'inlet'))]))
    check('Isolated COPV leaves the tank at its own pressure [psia]',
          isolated['blocks']['tank']['calculated']['calculated_pressure'], 50.0, 1e-12)
    through_closed_coupling = solve_model(small_model(
        pressurant_blocks + [('bottle', 'copv', 0, {'ullage_pressure_setpoint': 300.0}),
                             ('coupling', 'fluid_coupling', 1, {})],
        feed + [(('bottle', 'outlet'), ('coupling', 'inlet')), (('coupling', 'outlet'), ('tank', 'inlet'))]))
    check('Closed coupling blocks COPV backfill [psia]',
          through_closed_coupling['blocks']['tank']['calculated']['calculated_pressure'], 50.0, 1e-12)
    # open coupling in the feed: flow follows the two equivalent orifices in series [SUTTON]
    coupling_feed = small_model([('tank', 'tank', 0, {'nominal_pressure': 300.0}),
                                 ('coupling', 'fluid_coupling', 0, {'orifice_diameter': 0.02}),
                                 ('engine', 'thruster', 0, {'inlet_tubing_size': DEFAULT_TUBE})],
                                [(('tank', 'outlet'), ('coupling', 'inlet')),
                                 (('coupling', 'outlet'), ('engine', 'inlet'))])
    for line in coupling_feed['lines']:
        line['properties']['length'] = 0.001
    coupling_results = solve_model(coupling_feed)
    coupling_area = 0.6 * math.pi * (0.02 * METERS_PER_INCH) ** 2 / 4
    thruster_area = 0.6 * math.pi * (0.0143 * METERS_PER_INCH) ** 2 / 4
    series_area = 1 / math.sqrt(1 / coupling_area ** 2 + 1 / thruster_area ** 2)
    check('Open coupling: series equivalent-orifice flow [kg/s]', coupling_results['blocks']['engine']['mass_flow'],
          series_area * math.sqrt(2 * DEFAULT_FLUID['density'] * 300 * PASCALS_PER_PSI), 1e-3)
    dead_end = solve_model(small_model([('tank', 'tank', 0, {'nominal_pressure': 300.0}),
                                        ('service', 'fill_drain_valve', 0, {}),
                                        ('engine', 'thruster', 0, {'inlet_tubing_size': DEFAULT_TUBE})],
                                       [(('tank', 'outlet'), ('service', 'inlet')),
                                        (('service', 'outlet'), ('engine', 'inlet'))]))
    check('Fill & drain valve passes no flow [kg/s]', dead_end['blocks']['engine']['mass_flow'] + 1.0, 1.0, 1e-12)
    check('Fill & drain valve reads the connected pressure [psia]',
          dead_end['blocks']['service']['calculated']['calculated_inlet_pressure'], 300.0, 1e-9)

    print("Copy names")
    taken = {'LV-1', 'LV-2', 'Main valve', 'TK-09'}
    naming_cases = [('LV-1', 'LV-3'), ('LV-2', 'LV-3'), ('Main valve', 'Main valve 1'), ('TK-09', 'TK-10'),
                    ('Valve 2 left', 'Valve 3 left'), ('Pump-', 'Pump-1')]
    naming_ok = all(incremented_name(original, taken) == expected for original, expected in naming_cases)
    print(f"  [{'PASS' if naming_ok else 'FAIL'}] copy names increment the last number "
          f"({', '.join(f'{original} -> {incremented_name(original, taken)}' for original, _ in naming_cases)})")
    all_passed &= naming_ok
    print('ALL PASSED' if all_passed else 'SOME CHECKS FAILED')
    return 0 if all_passed else 1


if __name__ == '__main__' and '--selftest' in sys.argv:
    sys.exit(run_self_test())


# =============================================================================
# GUI (PySide6)
# =============================================================================
try:
    from PySide6.QtCore import QMimeData, QPointF, QRectF, Qt, QTimer, Signal
    from PySide6.QtGui import (QAction, QActionGroup, QBrush, QColor, QDrag, QFont,
                               QFontMetricsF, QImage, QKeySequence, QLinearGradient, QPainter,
                               QPainterPath, QPainterPathStroker, QPalette, QPen, QPixmap,
                               QPolygonF)
    from PySide6.QtWidgets import (QAbstractItemView, QAbstractSpinBox, QApplication, QComboBox, QDialog,
                                   QDialogButtonBox, QDockWidget, QFileDialog, QFormLayout,
                                   QFrame, QGraphicsItem, QGraphicsPathItem, QGraphicsScene,
                                   QGraphicsTextItem, QGraphicsView, QGraphicsEllipseItem,
                                   QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit,
                                   QListWidget, QListWidgetItem, QMainWindow, QMenu,
                                   QMessageBox, QProgressDialog, QPushButton, QScrollArea,
                                   QSpinBox, QSplitter, QTabWidget, QTableWidget,
                                   QTableWidgetItem, QTextEdit, QPlainTextEdit, QToolTip, QVBoxLayout,
                                   QWidget)
except ImportError:
    if __name__ == '__main__':
        sys.stderr.write("PySide6 is required for the GUI:\n    pip install PySide6 numpy\n"
                         "(run with --selftest to check the solver without the GUI)\n")
        sys.exit(1)
    raise

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
BLOCK_DIRECTORY = os.path.join(SCRIPT_DIRECTORY, 'Component Blocks')
BLOCK_MIME_TYPE = 'application/x-pylon-block'

# -----------------------------------------------------------------------------
# Theme: dark panels in gold bezels with a red inner trim and cyan accents
# -----------------------------------------------------------------------------
THEME = {
    'background_deep': '#04060c', 'background_middle': '#0a1020', 'background_panel': '#0b1424',
    'background_panel_top': '#14243c', 'background_field': '#060b14', 'background_alternate_row': '#0a1322',
    'gold_highlight': '#ffe9a8', 'gold': '#d9a441', 'gold_middle': '#a8741f', 'gold_dark': '#5a3a0c',
    'trim': '#b8361f',
    'cyan': '#5fd4ff', 'cyan_soft': '#8fd3ff', 'cyan_dim': '#2a6f8f',
    'text': '#e6f4ff', 'muted': '#8fa9c0', 'heading': '#f2c14e',
    'go': '#7dff8a', 'warning': '#ff6b5a', 'success': '#66e07a',
    'row_warning': '#3a1414', 'row_good': '#123a24', 'row_dimmed': '#0d1626',
    'off': '#4d5b6b',
}
COLOR_LINE = QColor(THEME['gold'])
COLOR_LINE_SELECTED = QColor(THEME['cyan'])
COLOR_WARNING = QColor(THEME['warning'])
COLOR_SELECTED = QColor(THEME['cyan'])
COLOR_TEXT = QColor(THEME['text'])
COLOR_MUTED = QColor(THEME['muted'])
COLOR_HEADING = QColor(THEME['heading'])
COLOR_SYMBOL = QColor('#7fe0ff')
FONT_FAMILIES = ['Bahnschrift', 'Segoe UI', 'Helvetica Neue', 'DejaVu Sans', 'Arial']


def draw_bezel(painter, rectangle, radius=12.0, rim=6.0, fill=True, glow=None, title=None,
               dimmed=False, gem=None, gem_side='top'):
    """Gold bevelled frame: metallic rim, red inner trim, dark panel, corner fins.

    Used for every panel, dialog and component block so the whole UI shares
    one frame language.  `glow` (QColor) adds a soft halo (selection/warning).
    """
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    outer_rectangle = QRectF(rectangle)
    if glow is not None:
        for ring_index in range(5, 0, -1):
            halo_color = QColor(glow)
            halo_color.setAlpha(int(150 / (ring_index + 1)))
            painter.setPen(QPen(halo_color, 2.0))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRoundedRect(outer_rectangle.adjusted(-ring_index, -ring_index, ring_index, ring_index),
                                    radius + ring_index, radius + ring_index)
    outer_path = QPainterPath()
    outer_path.addRoundedRect(outer_rectangle, radius, radius)
    inner_rectangle = outer_rectangle.adjusted(rim, rim, -rim, -rim)
    inner_radius = max(radius - rim * 0.6, 2.0)
    inner_path = QPainterPath()
    inner_path.addRoundedRect(inner_rectangle, inner_radius, inner_radius)
    if fill:
        panel_gradient = QLinearGradient(inner_rectangle.topLeft(), inner_rectangle.bottomLeft())
        if dimmed:
            panel_gradient.setColorAt(0.0, QColor('#171b22'))
            panel_gradient.setColorAt(1.0, QColor('#0a0c10'))
        else:
            panel_gradient.setColorAt(0.0, QColor(THEME['background_panel_top']))
            panel_gradient.setColorAt(0.45, QColor(THEME['background_panel']))
            panel_gradient.setColorAt(1.0, QColor('#060a12'))
        painter.fillPath(inner_path, QBrush(panel_gradient))
    # corner fins (drawn first so the rim overlaps their roots)
    fin_size = max(rim * 0.9, 3.0)
    fin_gradient = QLinearGradient(outer_rectangle.topLeft(), outer_rectangle.bottomRight())
    fin_gradient.setColorAt(0.0, QColor(THEME['gold_highlight']))
    fin_gradient.setColorAt(1.0, QColor(THEME['gold_middle']))
    painter.setPen(QPen(QColor(30, 18, 2, 220), 1.0))
    painter.setBrush(QBrush(fin_gradient))
    corner_inset = radius * 0.55
    for corner_x, corner_y, sign_x, sign_y in ((outer_rectangle.left(), outer_rectangle.top(), 1, 1),
                                               (outer_rectangle.right(), outer_rectangle.top(), -1, 1),
                                               (outer_rectangle.left(), outer_rectangle.bottom(), 1, -1),
                                               (outer_rectangle.right(), outer_rectangle.bottom(), -1, -1)):
        painter.drawPolygon(QPolygonF([
            QPointF(corner_x + sign_x * (corner_inset + fin_size * 1.6), corner_y + sign_y * corner_inset * 0.35),
            QPointF(corner_x - sign_x * fin_size * 0.55, corner_y - sign_y * fin_size * 0.55),
            QPointF(corner_x + sign_x * corner_inset * 0.35, corner_y + sign_y * (corner_inset + fin_size * 1.6))]))
    rim_path = outer_path.subtracted(inner_path)
    rim_gradient = QLinearGradient(outer_rectangle.topLeft(), outer_rectangle.bottomRight())
    if dimmed:
        for position, color in ((0, '#c9b07a'), (0.3, '#8a6a2c'), (0.6, '#b8975a'), (1, '#4a3510')):
            rim_gradient.setColorAt(position, QColor(color))
    else:
        for position, color in ((0, THEME['gold_highlight']), (0.18, THEME['gold']), (0.36, THEME['gold_dark']),
                                (0.55, '#f0cf7a'), (0.75, THEME['gold_middle']), (1, THEME['gold_dark'])):
            rim_gradient.setColorAt(position, QColor(color))
    painter.fillPath(rim_path, QBrush(rim_gradient))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.setPen(QPen(QColor(25, 14, 2, 230), 1.2))
    painter.drawPath(outer_path)
    painter.setPen(QPen(QColor(255, 240, 190, 140), 1.0))
    painter.drawRoundedRect(outer_rectangle.adjusted(1.4, 1.4, -1.4, -1.4), radius - 1.2, radius - 1.2)
    painter.setPen(QPen(QColor(THEME['trim']), 1.6))
    painter.drawPath(inner_path)
    painter.setPen(QPen(QColor(0, 0, 0, 170), 1.0))
    painter.drawRoundedRect(inner_rectangle.adjusted(1.5, 1.5, -1.5, -1.5),
                            max(inner_radius - 1.5, 1.0), max(inner_radius - 1.5, 1.0))
    if gem is not None:
        gem_size = max(rim * 0.75, 3.0)
        gem_x, gem_y = {'top': (outer_rectangle.center().x(), outer_rectangle.top() + rim / 2),
                        'bottom': (outer_rectangle.center().x(), outer_rectangle.bottom() - rim / 2),
                        'left': (outer_rectangle.left() + rim / 2, outer_rectangle.center().y()),
                        'right': (outer_rectangle.right() - rim / 2, outer_rectangle.center().y())}[gem_side]
        half_width, half_height = ((gem_size * 1.3, gem_size) if gem_side in ('top', 'bottom')
                                   else (gem_size, gem_size * 1.3))
        painter.setPen(QPen(QColor(THEME['gold_dark']), 1.0))
        painter.setBrush(QColor(gem))
        painter.drawPolygon(QPolygonF([QPointF(gem_x, gem_y - half_height), QPointF(gem_x + half_width, gem_y),
                                       QPointF(gem_x, gem_y + half_height), QPointF(gem_x - half_width, gem_y)]))
    if title:
        title_font = make_font(9.0, True)
        title_font.setLetterSpacing(QFont.SpacingType.PercentageSpacing, 112)
        painter.setFont(title_font)
        metrics = QFontMetricsF(title_font)
        text_x, text_y = inner_rectangle.left() + 10, inner_rectangle.top() + 4
        painter.setPen(QColor(0, 0, 0, 200))
        painter.drawText(QPointF(text_x + 1, text_y + metrics.ascent() + 1), title.upper())
        painter.setPen(QColor(THEME['heading']))
        painter.drawText(QPointF(text_x, text_y + metrics.ascent()), title.upper())
        rule_y = text_y + metrics.height() + 2
        rule_gradient = QLinearGradient(QPointF(text_x, rule_y), QPointF(inner_rectangle.right() - 10, rule_y))
        rule_gradient.setColorAt(0.0, QColor(THEME['gold']))
        rule_gradient.setColorAt(1.0, QColor(217, 164, 65, 0))
        painter.setPen(QPen(QBrush(rule_gradient), 1.2))
        painter.drawLine(QPointF(text_x, rule_y), QPointF(inner_rectangle.right() - 10, rule_y))
    painter.restore()


APPLICATION_STYLE_SHEET = """
* { font-family: "Bahnschrift", "Segoe UI", "Helvetica Neue", "DejaVu Sans", sans-serif; font-size: 9pt; }
QMainWindow { background: qradialgradient(cx:0.55, cy:0.35, radius:1.1, fx:0.55, fy:0.3,
              stop:0 #1b1546, stop:0.45 #0a0f1f, stop:1 #020308); }
QWidget { color: %(text)s; }
QDialog { background: %(background_middle)s; }
QLabel { background: transparent; }
QToolTip { background: %(background_panel)s; color: %(text)s; border: 1px solid %(gold)s; padding: 4px; }
QMenuBar { background: #070b14; color: %(heading)s; border-bottom: 2px solid %(gold_middle)s; padding: 2px; }
QMenuBar::item { background: transparent; padding: 4px 10px; }
QMenuBar::item:selected { background: #1a2a44; color: #ffffff; border-radius: 4px; }
QMenu { background: %(background_panel)s; border: 2px solid %(gold)s; color: %(text)s; padding: 4px; }
QMenu::item { padding: 5px 22px 5px 18px; border-radius: 4px; }
QMenu::item:selected { background: #1d3b5c; color: #ffffff; }
QMenu::item:disabled { color: %(off)s; }
QMenu::separator { height: 1px; background: %(gold_dark)s; margin: 4px 8px; }
QToolBar { background: #070b14; border: none; border-bottom: 2px solid %(gold_middle)s; spacing: 6px; padding: 5px; }
QToolBar::separator { background: %(gold_dark)s; width: 2px; margin: 4px 6px; }
QToolButton, QPushButton {
    color: %(text)s; padding: 4px 12px; min-height: 20px;
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #1a2c48, stop:1 #08101c);
    border: 2px solid %(gold)s; border-radius: 9px; }
QToolButton:hover, QPushButton:hover { border-color: %(gold_highlight)s; color: #ffffff;
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #22406a, stop:1 #0b1626); }
QToolButton:pressed, QPushButton:pressed { background: #050910; }
QToolButton:checked { color: %(cyan)s; border-color: %(cyan)s; background: #10304c; }
QPushButton:disabled, QToolButton:disabled { color: %(off)s; border-color: %(gold_dark)s; }
QPushButton[accent="go"] { color: %(go)s; font-weight: 600; }
QPushButton[accent="stop"] { color: %(warning)s; font-weight: 600; }
QLineEdit, QComboBox, QSpinBox {
    background: %(background_field)s; color: %(text)s; border: 1px solid %(gold_middle)s; border-radius: 5px;
    padding: 3px 6px; selection-background-color: #1d5f8a; }
QLineEdit:focus, QComboBox:focus, QSpinBox:focus { border-color: %(cyan)s; }
QLineEdit[readOnly="true"] { background: #0d1626; color: %(cyan_soft)s; border-color: #2a3a4c; }
QComboBox::drop-down { border: none; width: 18px; }
QComboBox QAbstractItemView { background: %(background_panel)s; color: %(text)s; border: 1px solid %(gold)s;
    selection-background-color: #1d3b5c; }
QTableWidget, QListWidget {
    background: %(background_field)s; alternate-background-color: %(background_alternate_row)s; color: %(text)s;
    gridline-color: #1b2a3d; border: 1px solid %(gold_dark)s; border-radius: 6px;
    selection-background-color: #1d3b5c; selection-color: #ffffff; }
QHeaderView::section { background: #0f1a2c; color: %(heading)s; border: none;
    border-bottom: 1px solid %(gold_middle)s; border-right: 1px solid #1b2a3d; padding: 4px 6px; font-weight: 600; }
QTableCornerButton::section { background: #0f1a2c; border: none; }
QDockWidget { color: %(heading)s; }
QDockWidget > QWidget { background: transparent; }
QTabBar::tab { background: %(background_panel)s; color: %(muted)s; border: 2px solid %(gold_dark)s; border-bottom: none;
    border-top-left-radius: 9px; border-top-right-radius: 9px; padding: 5px 16px; margin-right: 3px; }
QTabBar::tab:selected { color: %(heading)s; border-color: %(gold)s; background: #13223a; }
QTabBar::tab:hover { color: #ffffff; }
QScrollArea, QScrollArea > QWidget > QWidget { background: transparent; border: none; }
QScrollBar:vertical { background: %(background_field)s; width: 12px; margin: 0; border: 1px solid #2a1d06; border-radius: 6px; }
QScrollBar::handle:vertical { min-height: 26px; border-radius: 5px;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 %(gold_dark)s, stop:0.5 %(gold)s, stop:1 %(gold_dark)s); }
QScrollBar:horizontal { background: %(background_field)s; height: 12px; margin: 0; border: 1px solid #2a1d06; border-radius: 6px; }
QScrollBar::handle:horizontal { min-width: 26px; border-radius: 5px;
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 %(gold_dark)s, stop:0.5 %(gold)s, stop:1 %(gold_dark)s); }
QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: none; }
QStatusBar { background: #070b14; color: %(cyan_soft)s; border-top: 2px solid %(gold_middle)s; }
QSplitter::handle { background: transparent; }
QProgressBar { border: 1px solid %(gold_middle)s; border-radius: 5px; background: %(background_field)s; color: %(text)s; text-align: center; }
QProgressBar::chunk { background: %(cyan_dim)s; border-radius: 4px; }
QGraphicsView { background: %(background_deep)s; border: none; }
QRubberBand { border: 1px solid %(cyan)s; background: rgba(95, 212, 255, 40); }
QCheckBox, QRadioButton { color: %(text)s; }
QMessageBox QLabel, QInputDialog QLabel { color: %(text)s; }
""" % THEME


def apply_application_theme(application):
    application.setStyle('Fusion')
    palette = QPalette()
    for role, color in ((QPalette.ColorRole.Window, THEME['background_middle']),
                        (QPalette.ColorRole.WindowText, THEME['text']),
                        (QPalette.ColorRole.Base, THEME['background_field']),
                        (QPalette.ColorRole.AlternateBase, THEME['background_alternate_row']),
                        (QPalette.ColorRole.Text, THEME['text']), (QPalette.ColorRole.Button, '#13233b'),
                        (QPalette.ColorRole.ButtonText, THEME['text']), (QPalette.ColorRole.Highlight, '#1d5f8a'),
                        (QPalette.ColorRole.HighlightedText, '#ffffff'),
                        (QPalette.ColorRole.ToolTipBase, THEME['background_panel']),
                        (QPalette.ColorRole.ToolTipText, THEME['text']),
                        (QPalette.ColorRole.PlaceholderText, '#5f7488'),
                        (QPalette.ColorRole.Link, THEME['cyan']), (QPalette.ColorRole.BrightText, '#ffffff')):
        palette.setColor(role, QColor(color))
    application.setPalette(palette)
    application.setStyleSheet(APPLICATION_STYLE_SHEET)


def accent_buttons(button_box):
    """Green text for accept/apply buttons and red for cancel/close, like 'Start'/'Cancel'."""
    for button in button_box.buttons():
        role = button_box.buttonRole(button)
        if role in (QDialogButtonBox.ButtonRole.AcceptRole, QDialogButtonBox.ButtonRole.ApplyRole,
                    QDialogButtonBox.ButtonRole.YesRole):
            button.setProperty('accent', 'go')
        elif role in (QDialogButtonBox.ButtonRole.RejectRole, QDialogButtonBox.ButtonRole.NoRole):
            button.setProperty('accent', 'stop')
    return button_box


def to_qpoint(point):
    return QPointF(float(point[0]), float(point[1]))


def to_tuple(qpoint):
    return (qpoint.x(), qpoint.y())


def points_close(first, second, tolerance=0.6):
    return (first is not None and second is not None
            and abs(first[0] - second[0]) < tolerance and abs(first[1] - second[1]) < tolerance)


def distance_between(first, second):
    return math.hypot(first[0] - second[0], first[1] - second[1])


def make_font(point_size=8.0, bold=False):
    font = QFont()
    font.setFamilies(FONT_FAMILIES)
    font.setPointSizeF(point_size)
    font.setBold(bold)
    return font


def muted_html(text):
    return f"<span style='color:{THEME['muted']}'>{text}</span>"


def warning_html(text):
    return f"<span style='color:{THEME['warning']}'>{text}</span>"


class BezelFrame(QFrame):
    """A panel wrapped in the gold bezel, with an optional engraved title."""

    def __init__(self, title=None, child=None, rim=7.0, radius=14.0, padding=8):
        super().__init__()
        self.title = title
        self.rim = rim
        self.radius = radius
        margin = int(rim) + 4 + padding
        self.content_layout = QVBoxLayout(self)
        self.content_layout.setContentsMargins(margin, margin + (20 if title else 0), margin, margin)
        self.content_layout.setSpacing(6)
        if child is not None:
            self.content_layout.addWidget(child)

    def paintEvent(self, event):
        painter = QPainter(self)
        draw_bezel(painter, QRectF(self.rect()).adjusted(4, 4, -4, -4), self.radius, self.rim, title=self.title)
        painter.end()


class BezelDialog(QDialog):
    """Dialog painted as a bezel-framed briefing panel titled with the window title."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setContentsMargins(12, 30, 12, 10)

    def paintEvent(self, event):
        painter = QPainter(self)
        background = QLinearGradient(QPointF(0, 0), QPointF(0, self.height()))
        background.setColorAt(0.0, QColor('#101a30'))
        background.setColorAt(1.0, QColor('#03050a'))
        painter.fillRect(self.rect(), QBrush(background))
        draw_bezel(painter, QRectF(self.rect()).adjusted(5, 5, -5, -5), 16, 6, fill=False,
                   title=self.windowTitle())
        painter.end()


# -----------------------------------------------------------------------------
# Symbols
# -----------------------------------------------------------------------------
def draw_pipe(painter, start, end, body_color=None):
    """A short brass tube, drawn like the plumbing lines (dark outline, gold body, highlight)."""
    body_color = body_color or COLOR_LINE
    painter.save()
    painter.setBrush(Qt.BrushStyle.NoBrush)
    for width, color in ((7.5, QColor(20, 12, 2, 235)), (4.2, body_color), (1.2, QColor(255, 244, 210, 150))):
        pipe_pen = QPen(color, width)
        pipe_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        painter.setPen(pipe_pen)
        painter.drawLine(start, end)
    painter.restore()


def draw_symbol(painter, block_type, state, color):
    """Draw a component symbol centred at the origin in the block's base
    (unrotated) frame.  Flow runs left to right (the tank: top to bottom)."""
    base_width, base_height = BLOCK_TYPES[block_type]['size']
    half_width = base_width / 2.0
    dimmed = state != 0 and block_type in ('flowmeter', 'transducer', 'pump', 'thruster',
                                           'orifice', 'tee', 'elbow', 'copv', 'fluid_coupling')
    stroke_color = QColor(THEME['off']) if dimmed else QColor(color)
    pen = QPen(stroke_color, 2.0)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    translucent_fill = QColor(stroke_color)
    translucent_fill.setAlpha(70)
    interior_fill = QColor(THEME['background_field'])

    def flow_line(start_x, end_x):
        painter.drawLine(QPointF(start_x, 0), QPointF(end_x, 0))

    def draw_letters(rectangle, letters, point_size=8.0):
        painter.save()
        painter.setFont(make_font(point_size, True))
        painter.drawText(rectangle, Qt.AlignmentFlag.AlignCenter, letters)
        painter.restore()

    if block_type == 'tank':
        capsule = QRectF(-88, -160, 176, 280)              # 4x the earlier 44 x 70 capsule
        capsule_path = QPainterPath()
        capsule_path.addRoundedRect(capsule, 80, 80)
        fill_level = {0: 0.6, 1: 0.85, 2: 0.3}.get(state, 0.6)
        painter.save()
        painter.setClipPath(capsule_path)
        painter.fillRect(QRectF(capsule.left(), capsule.bottom() - capsule.height() * fill_level,
                                capsule.width(), capsule.height() * fill_level), translucent_fill)
        painter.restore()
        # inlet tube on top, outlet tube on the bottom (brass, like the plumbing lines)
        draw_pipe(painter, QPointF(0, -base_height / 2), QPointF(0, capsule.top() + 2))
        draw_pipe(painter, QPointF(0, capsule.bottom() - 2), QPointF(0, base_height / 2))
        tank_pen = QPen(stroke_color, 3.0)
        tank_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(tank_pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPath(capsule_path)
        painter.drawLine(QPointF(-14, capsule.bottom()), QPointF(14, capsule.bottom()))   # outlet boss
        painter.setPen(pen)
    elif block_type == 'copv':
        # composite overwrapped pressure vessel: capsule with cross-hatched overwrap, outlet tube below
        bottle = QRectF(-30, -55, 60, 110)
        bottle_path = QPainterPath()
        bottle_path.addRoundedRect(bottle, 28, 28)
        painter.save()
        painter.setClipPath(bottle_path)
        painter.fillRect(bottle, translucent_fill if state == 0 else interior_fill)
        painter.restore()
        hatch_color = QColor(stroke_color)
        hatch_color.setAlpha(110)
        painter.setPen(QPen(hatch_color, 1.0))
        corner_radius, half_width_inner, half_height_inner = 28.0, 30.0 - 28.0, 55.0 - 28.0

        def inside_bottle(x, y):
            dx, dy = max(abs(x) - half_width_inner, 0.0), max(abs(y) - half_height_inner, 0.0)
            return dx * dx + dy * dy <= (corner_radius - 1.5) ** 2

        for offset in range(-84, 85, 12):                        # +/-45 degree composite wrap
            for slope in (1, -1):
                segment_start = None
                for step in range(0, 121):                       # walk the line, keep the part inside
                    y = -60.0 + step
                    x = offset + slope * y
                    if inside_bottle(x, y):
                        if segment_start is None:
                            segment_start = (x, y)
                        segment_end = (x, y)
                    elif segment_start is not None:
                        painter.drawLine(QPointF(*segment_start), QPointF(*segment_end))
                        segment_start = None
                if segment_start is not None:
                    painter.drawLine(QPointF(*segment_start), QPointF(*segment_end))
        draw_pipe(painter, QPointF(0, bottle.bottom() - 2), QPointF(0, base_height / 2))
        painter.setPen(QPen(stroke_color, 2.5))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPath(bottle_path)
        painter.setPen(pen)
    elif block_type == 'fluid_coupling':
        # coupling halves meeting in the middle (open) or parted with sealed faces (closed)
        if state == 0:
            flow_line(-half_width, -20)
            flow_line(20, half_width)
            painter.setBrush(interior_fill)
            painter.drawRect(QRectF(-20, -9, 18, 18))
            painter.drawRect(QRectF(2, -9, 18, 18))
            painter.drawLine(QPointF(-2, 0), QPointF(2, 0))
            painter.setBrush(translucent_fill)
            painter.drawPolygon(QPolygonF([QPointF(-8, -5), QPointF(-2, 0), QPointF(-8, 5)]))
            painter.drawPolygon(QPolygonF([QPointF(8, -5), QPointF(2, 0), QPointF(8, 5)]))
        else:
            flow_line(-half_width, -24)
            flow_line(24, half_width)
            painter.setBrush(interior_fill)
            painter.drawRect(QRectF(-24, -9, 16, 18))
            painter.drawRect(QRectF(8, -9, 16, 18))
            sealed_pen = QPen(stroke_color, 3.5)
            painter.setPen(sealed_pen)
            painter.drawLine(QPointF(-8, -9), QPointF(-8, 9))
            painter.drawLine(QPointF(8, -9), QPointF(8, 9))
            painter.setPen(pen)
    elif block_type == 'fill_drain_valve':
        # normally closed hand valve (filled bow-tie) with a T handle: symbolic dead end
        flow_line(-half_width, -22)
        flow_line(22, half_width)
        painter.setBrush(stroke_color)
        painter.drawPolygon(QPolygonF([QPointF(-22, -13), QPointF(-22, 13), QPointF(0, 0)]))
        painter.drawPolygon(QPolygonF([QPointF(22, -13), QPointF(22, 13), QPointF(0, 0)]))
        painter.drawLine(QPointF(0, 0), QPointF(0, -17))
        painter.drawLine(QPointF(-9, -17), QPointF(9, -17))
    elif block_type == 'latch_valve':
        flow_line(-half_width, -22)
        flow_line(22, half_width)
        painter.setBrush(stroke_color if state == 1 else interior_fill)
        painter.drawPolygon(QPolygonF([QPointF(-22, -13), QPointF(-22, 13), QPointF(0, 0)]))
        painter.drawPolygon(QPolygonF([QPointF(22, -13), QPointF(22, 13), QPointF(0, 0)]))
        painter.setBrush(interior_fill)
        painter.drawLine(QPointF(0, 0), QPointF(0, -16))
        painter.drawRect(QRectF(-9, -25, 18, 9))
    elif block_type == 'filter':
        flow_line(-half_width, -20)
        flow_line(20, half_width)
        painter.setBrush(QBrush(stroke_color, Qt.BrushStyle.Dense5Pattern) if state == 1 else interior_fill)
        painter.drawPolygon(QPolygonF([QPointF(-20, 0), QPointF(0, -19), QPointF(20, 0), QPointF(0, 19)]))
        painter.setPen(QPen(stroke_color, 1.5, Qt.PenStyle.DashLine))
        painter.drawLine(QPointF(0, -19), QPointF(0, 19))
    elif block_type == 'flowmeter':
        flow_line(-half_width, -17)
        flow_line(17, half_width)
        painter.setBrush(interior_fill)
        painter.drawEllipse(QPointF(0, 0), 17, 17)
        draw_letters(QRectF(-17, -17, 34, 34), 'FM')
        if state == 1:
            painter.setPen(QPen(COLOR_WARNING, 2.5))
            painter.drawLine(QPointF(-12, -12), QPointF(12, 12))
            painter.drawLine(QPointF(-12, 12), QPointF(12, -12))
    elif block_type == 'transducer':
        painter.drawLine(QPointF(0, 9), QPointF(0, base_height / 2))
        painter.setBrush(interior_fill)
        painter.drawEllipse(QPointF(0, -8), 17, 17)
        draw_letters(QRectF(-17, -25, 34, 34), 'PT')
    elif block_type == 'pump':
        # gear pump symbol: casing circle, two touching gear circles, and a V whose point
        # sits on the casing at the outlet - the V points in the flow direction (left to right)
        casing_radius = 21.0
        flow_line(-half_width, -casing_radius)
        flow_line(casing_radius, half_width)
        painter.setBrush(interior_fill)
        painter.drawEllipse(QPointF(0, 0), casing_radius, casing_radius)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPolyline(QPolygonF([QPointF(0, -casing_radius), QPointF(casing_radius, 0),
                                        QPointF(0, casing_radius)]))
        gear_radius = 0.30 * casing_radius
        gear_center_x = -0.20 * casing_radius          # gears sit back from the V, toward the inlet
        painter.setBrush(translucent_fill if state == 0 else interior_fill)
        for gear_center_y in (-gear_radius, gear_radius):
            painter.drawEllipse(QPointF(gear_center_x, gear_center_y), gear_radius, gear_radius)
    elif block_type == 'thruster':
        flow_line(-half_width, -34)
        painter.setBrush(interior_fill)
        painter.drawRect(QRectF(-34, -11, 26, 22))
        painter.drawPolygon(QPolygonF([QPointF(-8, -6), QPointF(26, -21), QPointF(26, 21), QPointF(-8, 6)]))
        if state == 0:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(235, 104, 52, 150))
            painter.drawPolygon(QPolygonF([QPointF(27, -14), QPointF(52, 0), QPointF(27, 14)]))
    elif block_type == 'orifice':
        flow_line(-half_width, -4)
        flow_line(4, half_width)
        painter.setPen(QPen(stroke_color, 4.0))
        if state == 0:
            painter.drawLine(QPointF(0, -20), QPointF(0, -4))
            painter.drawLine(QPointF(0, 4), QPointF(0, 20))
        else:
            painter.drawLine(QPointF(0, -20), QPointF(0, 20))
    elif block_type == 'tee':
        painter.setPen(QPen(stroke_color, 5.0))
        flow_line(-half_width, half_width)
        painter.drawLine(QPointF(0, 0), QPointF(0, base_height / 2))
        painter.setBrush(stroke_color)
        painter.drawEllipse(QPointF(0, 0), 4, 4)
    elif block_type == 'cross':
        painter.setPen(QPen(stroke_color, 5.0))
        flow_line(-half_width, half_width)
        painter.drawLine(QPointF(0, -base_height / 2), QPointF(0, base_height / 2))
        painter.setBrush(stroke_color)
        painter.drawEllipse(QPointF(0, 0), 4, 4)
    elif block_type == 'elbow':
        painter.setPen(QPen(stroke_color, 5.0))
        elbow_path = QPainterPath(QPointF(-half_width, 0))
        elbow_path.lineTo(-14, 0)
        elbow_path.quadTo(QPointF(0, 0), QPointF(0, 14))
        elbow_path.lineTo(0, base_height / 2)
        painter.drawPath(elbow_path)


def symbol_pixmap(block_type, state=0, width=76, height=50):
    pixmap = QPixmap(width, height)
    pixmap.fill(QColor(0, 0, 0, 0))
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    base_width, base_height = BLOCK_TYPES[block_type]['size']
    scale = min(width / base_width, height / base_height) * 0.92
    painter.translate(width / 2, height / 2)
    painter.scale(scale, scale)
    draw_symbol(painter, block_type, state, COLOR_SYMBOL)
    painter.end()
    return pixmap


def assembly_pixmap(template, width=76, height=50):
    """Palette icon for a saved assembly: up to three of its parts side by side."""
    pixmap = QPixmap(width, height)
    pixmap.fill(QColor(0, 0, 0, 0))
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    part_types = [block.get('type') for block in template.get('blocks', []) if block.get('type') in BLOCK_TYPES][:3]
    slot_width = width / max(len(part_types), 1)
    for slot_index, block_type in enumerate(part_types):
        base_width, base_height = BLOCK_TYPES[block_type]['size']
        scale = min(slot_width / base_width, height / base_height) * 0.9
        painter.save()
        painter.translate(slot_width * (slot_index + 0.5), height / 2)
        painter.scale(scale, scale)
        draw_symbol(painter, block_type, 0, COLOR_SYMBOL)
        painter.restore()
    painter.setPen(QPen(QColor(THEME['gold_middle']), 1.0, Qt.PenStyle.DashLine))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRoundedRect(QRectF(1, 1, width - 2, height - 2), 6, 6)
    painter.end()
    return pixmap


# -----------------------------------------------------------------------------
# Small graphics items
# -----------------------------------------------------------------------------
class LabelItem(QGraphicsItem):
    """Text, optionally on a gold-edged plate; never takes mouse input."""

    def __init__(self, parent=None, point_size=7.5, bold=False, color=COLOR_TEXT, plate=True):
        super().__init__(parent)
        self.plate = plate
        self.text = ''
        self.font = make_font(point_size, bold)
        self.color = QColor(color)
        self.rectangle = QRectF()
        self.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        self.setZValue(10)

    def set_text(self, text):
        self.prepareGeometryChange()
        self.text = text
        if text:
            metrics = QFontMetricsF(self.font)
            rows = text.split('\n')
            width = max(metrics.horizontalAdvance(row) for row in rows) + 12
            height = metrics.height() * len(rows) + 4
            self.rectangle = QRectF(-width / 2, -height / 2, width, height)
        else:
            self.rectangle = QRectF()
        self.update()

    def boundingRect(self):
        return self.rectangle

    def paint(self, painter, option, widget=None):
        if not self.text:
            return
        if self.plate:
            painter.setPen(QPen(QColor(217, 164, 65, 170), 1.0))
            painter.setBrush(QColor(6, 11, 20, 230))
            painter.drawRoundedRect(self.rectangle, 4, 4)
        else:
            painter.setFont(self.font)
            painter.setPen(QColor(0, 0, 0, 220))
            painter.drawText(self.rectangle.translated(1, 1), Qt.AlignmentFlag.AlignCenter, self.text)
        painter.setPen(self.color)
        painter.setFont(self.font)
        painter.drawText(self.rectangle, Qt.AlignmentFlag.AlignCenter, self.text)


class PortItem(QGraphicsEllipseItem):
    RADIUS = 5.0

    def __init__(self, block, name):
        super().__init__(-self.RADIUS, -self.RADIUS, 2 * self.RADIUS, 2 * self.RADIUS, block)
        self.block = block
        self.name = name
        self.direction = (1, 0)
        self.connected = False
        self.setZValue(3)
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.CursorShape.CrossCursor)
        self.setToolTip(f"{name} port — drag to another port or line to connect")
        self._apply_style(False)

    def _apply_style(self, hovered):
        self.setPen(QPen(QColor(THEME['gold_highlight']), 1.6))
        if hovered:
            self.setBrush(QColor(THEME['cyan']))
        else:
            self.setBrush(QColor(THEME['gold']) if self.connected else QColor(THEME['background_deep']))

    def set_connected(self, connected):
        self.connected = connected
        self._apply_style(False)

    def hoverEnterEvent(self, event):
        self._apply_style(True)

    def hoverLeaveEvent(self, event):
        self._apply_style(False)

    def mousePressEvent(self, event):
        scene = self.scene()
        if (event.button() == Qt.MouseButton.LeftButton and isinstance(scene, SchematicScene)
                and scene.mode == 'select'):
            scene.begin_connection(self)
            event.accept()
        else:
            event.ignore()

    def mouseMoveEvent(self, event):
        scene = self.scene()
        if isinstance(scene, SchematicScene) and scene.connection:
            scene.update_connection(event.scenePos())

    def mouseReleaseEvent(self, event):
        scene = self.scene()
        if isinstance(scene, SchematicScene) and scene.connection:
            scene.end_connection(event.scenePos())


class SpeedEditorItem(QGraphicsTextItem):
    """In-block editor for the pump speed (rpm), shown while the pump is ON."""

    def __init__(self, block):
        super().__init__(block)
        self.block = block
        self.setFont(make_font(8.0, True))
        self.setDefaultTextColor(COLOR_TEXT)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        self.setToolTip("Pump speed [rpm] — click to edit, Enter to apply. The pump follows its "
                        "imported curve for this speed (interpolated between stored speeds).")
        self.setZValue(4)

    def synchronize(self):
        self.setPlainText(f"{float(self.block.properties.get('revolutions_per_minute', 0)):.0f}")

    def keyPressEvent(self, event):
        if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self.clearFocus()
            return
        if event.key() == Qt.Key.Key_Escape:
            self.synchronize()
            self.clearFocus()
            return
        super().keyPressEvent(event)

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self.commit()

    def commit(self):
        try:
            value = float(self.toPlainText().strip())
        except ValueError:
            value = None
        if (value is not None and value >= 0
                and abs(value - float(self.block.properties.get('revolutions_per_minute', 0))) > 1e-9):
            self.block.properties['revolutions_per_minute'] = value
            scene = self.scene()
            if isinstance(scene, SchematicScene):
                scene.notify()
        self.synchronize()
        self.block.update()

    def paint(self, painter, option, widget=None):
        painter.setPen(QPen(QColor(THEME['cyan']), 1.0))
        painter.setBrush(QColor(THEME['background_field']))
        painter.drawRoundedRect(self.boundingRect().adjusted(1, 2, -1, -2), 3, 3)
        super().paint(painter, option, widget)


class AssemblyOutlineItem(QGraphicsItem):
    """Faint dashed outline and name tag around the parts of one assembly; a padlock shows
    whether the assembly is locked (moves and edits as one block) or unlocked."""

    def __init__(self, assembly_identifier):
        super().__init__()
        self.assembly_identifier = assembly_identifier
        self.rectangle = QRectF()
        self.name = ''
        self.locked = True
        self.setZValue(-3)
        self.setAcceptedMouseButtons(Qt.MouseButton.NoButton)

    def set_geometry(self, rectangle, name, locked):
        self.prepareGeometryChange()
        self.rectangle, self.name, self.locked = rectangle, name, locked
        self.update()

    def boundingRect(self):
        return self.rectangle.adjusted(-4, -24, 4, 4)

    def paint(self, painter, option, widget=None):
        if self.rectangle.width() <= 0 or self.rectangle.height() <= 0:
            return
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        color = QColor(THEME['gold_middle'] if self.locked else THEME['cyan_dim'])
        painter.setPen(QPen(color, 1.4, Qt.PenStyle.DashLine))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRoundedRect(self.rectangle, 12, 12)
        # padlock: body plus shackle (shackle lifted when unlocked)
        tag_x, tag_y = self.rectangle.left() + 6, self.rectangle.top() - 18
        painter.setPen(QPen(color, 1.4))
        painter.setBrush(color)
        painter.drawRect(QRectF(tag_x, tag_y + 6, 9, 7))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        lift = 0 if self.locked else 3
        painter.drawPolyline(QPolygonF([QPointF(tag_x + 2, tag_y + 6), QPointF(tag_x + 2, tag_y + 2 - lift),
                                        QPointF(tag_x + 7, tag_y + 2 - lift),
                                        QPointF(tag_x + 7, tag_y + 6 - lift)]))
        painter.setFont(make_font(8.0, True))
        painter.setPen(color)
        painter.drawText(QRectF(tag_x + 14, tag_y, 400, 16), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                         f"{self.name} · {'locked' if self.locked else 'unlocked'}")


# -----------------------------------------------------------------------------
# Flow / pressure data label visibility (shared by blocks and lines)
# -----------------------------------------------------------------------------
class DataLabelVisibility:
    """Flow / pressure labels follow the scene-wide switch (scene.show_labels) unless
    this element has its own override (True / False, set from its right-click menu).
    While the mouse is over a block or line its hidden labels are shown temporarily."""
    data_label_override = None
    _hover_reveal = False

    def data_labels_shown(self):
        """Visibility ignoring the hover reveal: the override, else the scene-wide switch."""
        if self.data_label_override is not None:
            return bool(self.data_label_override)
        return bool(getattr(self.scene(), 'show_labels', True))

    def set_data_label_override(self, visible):
        global_visible = bool(getattr(self.scene(), 'show_labels', True))
        self.data_label_override = None if bool(visible) == global_visible else bool(visible)
        self.apply_data_label_visibility()

    def set_hover_reveal(self, revealed):
        if revealed != self._hover_reveal:
            self._hover_reveal = revealed
            self.apply_data_label_visibility()

    def data_labels_visible_now(self):
        return self.data_labels_shown() or self._hover_reveal

    def apply_data_label_visibility(self):
        raise NotImplementedError


# -----------------------------------------------------------------------------
# Blocks
# -----------------------------------------------------------------------------
class BlockItem(DataLabelVisibility, QGraphicsItem):
    def __init__(self, data):
        super().__init__()
        self.identifier = data.get('identifier') or new_identifier()
        self.block_type = data['type']
        self.type_information = BLOCK_TYPES[self.block_type]
        self.name = data.get('name', '')
        self.properties = default_block_properties(self.block_type)
        self.properties.update(copy.deepcopy(data.get('properties') or {}))
        remove_obsolete_properties(self.block_type, self.properties)
        update_linked_flow_coefficients(self.block_type, self.properties)
        self.state = normalized_state(self.block_type, data.get('state', 0))
        self.pump_curves = copy.deepcopy(data.get('pump_curves') or {})
        self.orientation = int(data.get('orientation', 0)) % 360
        self.data_label_override = data.get('data_label_visibility')
        self._hover_reveal = False
        self.assembly_identifier = data.get('assembly')
        self._last_position = None
        self.result = None
        self.warnings = []
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setAcceptHoverEvents(True)
        self.setZValue(1)
        self.ports = {port_name: PortItem(self, port_name) for port_name in self.type_information['ports']}
        self.result_label = LabelItem(self, point_size=7.5)
        self.name_label = LabelItem(self, point_size=8.5, bold=True, color=COLOR_HEADING, plate=False)
        self.speed_editor = SpeedEditorItem(self) if self.block_type == 'pump' else None
        self._press_position = None
        self.setPos(float(data.get('x', 0)), float(data.get('y', 0)))
        self.layout_ports()
        self.refresh()

    # geometry -------------------------------------------------------------
    def size(self):
        return block_size(self.block_type, self.orientation)

    def body_rectangle(self):
        width, height = self.size()
        return QRectF(-width / 2, -height / 2, width, height)

    def scene_rectangle(self):
        return block_rectangle({'type': self.block_type, 'x': self.pos().x(), 'y': self.pos().y(),
                                'orientation': self.orientation})

    def layout_ports(self):
        self.prepareGeometryChange()
        for port_name, port in self.ports.items():
            (position_x, position_y), direction = port_local_position(self.block_type, port_name, self.orientation)
            port.setPos(position_x, position_y)
            port.direction = direction
        self.place_side_labels()

    # ---- labels inside the block's (invisible) box, close to the symbol ------------
    LABEL_GAP = 3.0
    CAPTION_HEIGHT = 13.0
    SPEED_ROW_ON_WIDTH, SPEED_ROW_UNIT_WIDTH = 16.0, 21.0

    def symbol_extent(self):
        """How far the symbol reaches up, down, left and right of the centre, after rotation."""
        up, down, left, right = self.type_information.get('symbol_extent', (20, 20, 20, 20))
        reach = {}
        for distance, (direction_x, direction_y) in ((up, (0, -1)), (down, (0, 1)), (left, (-1, 0)), (right, (1, 0))):
            rotated_x, rotated_y = rotate_vector(direction_x, direction_y, self.orientation)
            reach[(int(rotated_x), int(rotated_y))] = float(distance)
        return reach[(0, -1)], reach[(0, 1)], reach[(-1, 0)], reach[(1, 0)]

    def _caption_row_size(self):
        """(width, height) of the row under the symbol: the pump speed box row, the state
        caption, or nothing."""
        if self.speed_editor is not None and self.state == 0:
            editor_rectangle = self.speed_editor.boundingRect()
            return (self.SPEED_ROW_ON_WIDTH + editor_rectangle.width() + self.SPEED_ROW_UNIT_WIDTH + 6,
                    max(editor_rectangle.height(), self.CAPTION_HEIGHT))
        if has_toggle(self.block_type) and self.block_type != 'tank':
            return 90.0, self.CAPTION_HEIGHT
        return 0.0, 0.0

    def label_layout(self):
        """Centres of the name and result labels, and the caption row (left, top, width,
        height, alignment), in block coordinates.  Labels sit just outside the symbol,
        inside the block's box, on sides without ports."""
        width, height = self.size()
        gap = self.LABEL_GAP
        name_rectangle, result_rectangle = self.name_label.boundingRect(), self.result_label.boundingRect()
        name_width, name_height = name_rectangle.width(), name_rectangle.height()
        result_width, result_height = result_rectangle.width(), result_rectangle.height()
        row_width, row_height = self._caption_row_size()
        up, down, left, right = self.symbol_extent()
        center_y = (down - up) / 2                 # vertical middle of the symbol
        used_directions = {tuple(port.direction) for port in self.ports.values()}
        free = {side for side, direction in (('top', (0, -1)), ('bottom', (0, 1)), ('left', (-1, 0)),
                                             ('right', (1, 0))) if direction not in used_directions}
        layout = {'caption': None}

        def left_band(label_width, y):        # right-aligned just left of the symbol
            return -(left + gap) - label_width / 2, y

        def right_band(label_width, y):       # left-aligned just right of the symbol
            return right + gap + label_width / 2, y

        if self.block_type == 'tank':
            # inside the capsule: name, state caption, results stacked down the middle
            layout['name'] = (0.0, -62.0)
            layout['result'] = (0.0, 34.0)
            layout['caption'] = (-80.0, -26.0, 160.0, 24.0, 'center')
        elif 'top' in free and 'bottom' in free:
            # name just above the symbol; caption (or pump speed row) and results just below
            layout['name'] = (0.0, -(up + gap + name_height / 2))
            row_top = down + gap
            if row_height:
                layout['caption'] = (-row_width / 2, row_top, row_width, row_height, 'center')
                row_top += row_height + gap
            layout['result'] = (0.0, row_top + result_height / 2)
        elif 'top' in free or 'bottom' in free:
            above = 'top' in free
            sign = -1 if above else 1
            reach = up if above else down
            layout['name'] = (0.0, sign * (reach + gap + name_height / 2))
            if 'right' in free:
                layout['result'] = right_band(result_width, center_y)
                if row_height and 'left' in free:
                    layout['caption'] = (-(left + gap) - row_width, center_y - row_height / 2,
                                         row_width, row_height, 'right')
            elif 'left' in free:
                layout['result'] = left_band(result_width, center_y)
            else:
                # e.g. a tee: results between the symbol and the name
                result_y = sign * (reach + gap + result_height / 2)
                layout['result'] = (0.0, result_y)
                layout['name'] = (0.0, result_y + sign * (result_height / 2 + gap + name_height / 2))
        elif 'left' in free and 'right' in free:
            # rotated in-line blocks: name and caption left of the symbol, results right of it
            name_y = center_y - (name_height / 2 + 1 if row_height else 0)
            layout['name'] = left_band(name_width, name_y)
            if row_height:
                layout['caption'] = (-(left + gap) - row_width, name_y + name_height / 2 + 1,
                                     row_width, row_height, 'right')
            layout['result'] = right_band(result_width, center_y)
        elif 'left' in free or 'right' in free:
            band = left_band if 'left' in free else right_band
            layout['name'] = band(name_width, center_y - name_height / 2 - 1)
            layout['result'] = band(result_width, center_y + result_height / 2 + 1)
        else:
            # ports on all four edges (cross): name upper-left, results upper-right
            layout['name'] = (-(2 * gap) - name_width / 2, -(2 * gap) - name_height / 2)
            layout['result'] = (2 * gap + result_width / 2, -(2 * gap) - result_height / 2)
        return layout

    def place_side_labels(self):
        """Put the name and result labels (and the pump speed box) where label_layout() says."""
        if self.name_label.text != self.name:
            self.name_label.set_text(self.name)
        layout = self.label_layout()
        self.name_label.setPos(*layout['name'])
        self.result_label.setPos(*layout['result'])
        if self.speed_editor is not None:
            caption = layout['caption']
            if caption is not None:
                row_left, row_top, _, row_height, _ = caption
                editor_rectangle = self.speed_editor.boundingRect()
                self.speed_editor.setPos(row_left + self.SPEED_ROW_ON_WIDTH + 3,
                                         row_top + (row_height - editor_rectangle.height()) / 2)
        self.update()

    def label_rectangles(self):
        """Scene rectangles of the block body and its name/result labels."""
        left, top, right, bottom = self.scene_rectangle()
        rectangles = [(left - 4, top - 4, right + 4, bottom + 4)]
        for label in (self.name_label, self.result_label):
            if label.isVisible() and label.text:
                center, rectangle = label.scenePos(), label.boundingRect()
                rectangles.append((center.x() + rectangle.left(), center.y() + rectangle.top(),
                                   center.x() + rectangle.right(), center.y() + rectangle.bottom()))
        return rectangles

    def port_scene_position(self, port_name):
        return to_tuple(self.ports[port_name].scenePos())

    def port_direction(self, port_name):
        return self.ports[port_name].direction

    def boundingRect(self):
        width, height = self.size()
        return QRectF(-width / 2 - 8, -height / 2 - 8, width + 16, height + 16)

    def toggle_zone(self):
        width, height = self.size()
        return QRectF(-width * 0.28, -height * 0.32, width * 0.56, height * 0.56)

    # state ------------------------------------------------------------------
    def state_label(self):
        label = self.type_information['states'][self.state]
        if self.block_type == 'tank':
            backfill_copv = (self.result or {}).get('backfill_copv')
            if backfill_copv:
                label = f"COPV backfill · {format_number((self.result or {}).get('calculated', {}).get('calculated_pressure', float('nan')), 4)} psia"
            else:
                label = (f"{label.split()[0]} · "
                         f"{format_number(tank_pressure_psia(self.properties, self.state), 4)} psia")
        if self.block_type == 'pump' and self.state == 0:
            label = f"ON · {format_number(float(self.properties.get('revolutions_per_minute', 0)), 5)} rpm"
        return label

    def set_state(self, state_index):
        if not has_toggle(self.block_type):
            return
        self.state = int(state_index) % len(self.type_information['states'])
        self.refresh()
        scene = self.scene()
        if isinstance(scene, SchematicScene):
            scene.notify()

    def refresh(self):
        if self.speed_editor is not None:
            self.speed_editor.setVisible(self.state == 0)
            if not self.speed_editor.hasFocus():
                self.speed_editor.synchronize()
            self.layout_ports()
        if self.name_label.text != self.name:
            self.place_side_labels()
        tooltip_lines = [f"<b>{self.name}</b> — {self.type_information['label']}"]
        if has_toggle(self.block_type):
            tooltip_lines += [f"State: {self.type_information['states'][self.state]}",
                              "Click the centre to toggle state · right-click for menu · "
                              "double-click for properties"]
        else:
            tooltip_lines.append("Always open (no toggle) · right-click for menu · double-click for properties")
        if self.warnings:
            tooltip_lines += [warning_html(f"⚠ {warning}") for warning in self.warnings]
        self.setToolTip('<br>'.join(tooltip_lines))
        self.update()

    def set_result(self, result):
        self.result = result
        self.warnings = (result or {}).get('warnings', [])
        if result is None:
            self.result_label.set_text('')
        elif self.block_type == 'transducer':
            reading = result.get('reading', float('nan'))
            self.result_label.set_text('' if math.isnan(reading) else f"{format_number(reading, 4)} psia")
        else:
            self.result_label.set_text(result.get('label', ''))
        self.apply_data_label_visibility()
        self.place_side_labels()
        self.refresh()

    def to_dictionary(self):
        data = {'identifier': self.identifier, 'type': self.block_type, 'name': self.name,
                'x': self.pos().x(), 'y': self.pos().y(), 'orientation': self.orientation,
                'state': self.state, 'properties': copy.deepcopy(self.properties),
                'pump_curves': copy.deepcopy(self.pump_curves)}
        if self.data_label_override is not None:
            data['data_label_visibility'] = self.data_label_override
        if self.assembly_identifier:
            data['assembly'] = self.assembly_identifier
        return data

    def apply_data_label_visibility(self):
        visible = self.data_labels_visible_now()
        self.result_label.setVisible(visible and bool(self.result_label.text))
        # a label revealed by hovering is drawn above neighbouring blocks
        self.setZValue(5 if self._hover_reveal and not self.data_labels_shown() else 1)

    # painting ---------------------------------------------------------------
    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        width, height = self.size()
        body = QRectF(-width / 2, -height / 2, width, height)
        active = self.state == 0 or self.block_type in ('tank', 'filter')
        # no frame around components: selection / warning show as a soft translucent halo
        halo = COLOR_WARNING if self.warnings else (COLOR_SELECTED if self.isSelected() else None)
        if halo is not None:
            halo_fill = QColor(halo)
            halo_fill.setAlpha(46)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(halo_fill)
            painter.drawRoundedRect(body.adjusted(-5, -5, 5, 5), 14, 14)
        painter.save()
        painter.rotate(self.orientation)
        draw_symbol(painter, self.block_type, self.state, COLOR_SYMBOL)
        painter.restore()
        painter.setFont(make_font(7.2))
        painter.setPen(QColor(THEME['cyan_soft']) if active else QColor(THEME['off']))
        caption = self.label_layout()['caption']
        if caption is not None:
            row_left, row_top, row_width, row_height, alignment = caption
            if self.speed_editor is not None and self.speed_editor.isVisible():
                # pump speed row:  ON [rpm box] rpm
                editor_width = self.speed_editor.boundingRect().width()
                painter.drawText(QRectF(row_left, row_top, self.SPEED_ROW_ON_WIDTH, row_height),
                                 Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, 'ON')
                painter.drawText(QRectF(row_left + self.SPEED_ROW_ON_WIDTH + 6 + editor_width, row_top,
                                        self.SPEED_ROW_UNIT_WIDTH, row_height),
                                 Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, 'rpm')
            elif self.block_type == 'tank':
                painter.setFont(make_font(9.0, True))
                painter.drawText(QRectF(row_left, row_top, row_width, row_height), Qt.AlignmentFlag.AlignCenter,
                                 self.state_label())
            elif has_toggle(self.block_type):
                horizontal = {'center': Qt.AlignmentFlag.AlignHCenter, 'right': Qt.AlignmentFlag.AlignRight,
                              'left': Qt.AlignmentFlag.AlignLeft}[alignment]
                painter.drawText(QRectF(row_left, row_top, row_width, row_height),
                                 horizontal | Qt.AlignmentFlag.AlignVCenter, self.state_label())
        if self.warnings:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(COLOR_WARNING)
            painter.drawEllipse(QPointF(width / 2 - 11, -height / 2 + 11), 6, 6)
            painter.setPen(QColor('#ffffff'))
            painter.setFont(make_font(7.5, True))
            painter.drawText(QRectF(width / 2 - 17, -height / 2 + 5, 12, 12), Qt.AlignmentFlag.AlignCenter, '!')

    # interaction ------------------------------------------------------------
    def hoverEnterEvent(self, event):
        self.set_hover_reveal(True)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        self.set_hover_reveal(False)
        super().hoverLeaveEvent(event)

    def hoverMoveEvent(self, event):
        if has_toggle(self.block_type) and self.toggle_zone().contains(event.pos()):
            self.setCursor(Qt.CursorShape.PointingHandCursor)
        else:
            self.setCursor(Qt.CursorShape.SizeAllCursor)

    def mousePressEvent(self, event):
        self._press_position = None
        if (event.button() == Qt.MouseButton.LeftButton and has_toggle(self.block_type)
                and self.toggle_zone().contains(event.pos())):
            self._press_position = to_tuple(event.scenePos())
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        scene = self.scene()
        if (self._press_position is not None
                and distance_between(self._press_position, to_tuple(event.scenePos())) < 3):
            self._press_position = None
            self.set_state(self.state + 1)
        self._press_position = None
        if isinstance(scene, SchematicScene):
            scene.on_block_move_finished()

    def mouseDoubleClickEvent(self, event):
        scene = self.scene()
        if isinstance(scene, SchematicScene) and not (has_toggle(self.block_type)
                                                       and self.toggle_zone().contains(event.pos())):
            scene.main_window.edit_block(self)
        else:
            super().mouseDoubleClickEvent(event)

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionChange:
            return QPointF(snap_to_grid(value.x()), snap_to_grid(value.y()))
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            position = to_tuple(self.pos())
            previous = self._last_position or position
            self._last_position = position
            scene = self.scene()
            if isinstance(scene, SchematicScene):
                scene.on_block_moving(self, (position[0] - previous[0], position[1] - previous[1]))
        if change == QGraphicsItem.GraphicsItemChange.ItemSelectedHasChanged:
            scene = self.scene()
            if isinstance(scene, SchematicScene):
                scene.on_item_selection_changed(self)
        return super().itemChange(change, value)


# -----------------------------------------------------------------------------
# Plumbing lines
# -----------------------------------------------------------------------------
class LineItem(DataLabelVisibility, QGraphicsPathItem):
    def __init__(self, data):
        super().__init__()
        self.identifier = data.get('identifier') or new_identifier()
        self.name = data.get('name', '')
        self.points = [(float(point[0]), float(point[1])) for point in (data.get('points') or [])]
        self.start = copy.deepcopy(data.get('start'))
        self.end = copy.deepcopy(data.get('end'))
        self.properties = default_line_properties()
        self.properties.update(copy.deepcopy(data.get('properties') or {}))
        self.auto_routed = bool(data.get('auto_routed', True))
        self.assembly_identifier = data.get('assembly')
        self.data_label_override = data.get('data_label_visibility')
        self._hover_reveal = False
        self.result = None
        self.labels = []
        self.arrows = []
        self._drag = None
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(-1)
        self.setAcceptHoverEvents(True)
        self._hovered = False

    def to_dictionary(self):
        data = {'identifier': self.identifier, 'name': self.name,
                'points': [[float(point[0]), float(point[1])] for point in self.points],
                'start': copy.deepcopy(self.start), 'end': copy.deepcopy(self.end),
                'properties': copy.deepcopy(self.properties), 'auto_routed': self.auto_routed}
        if self.data_label_override is not None:
            data['data_label_visibility'] = self.data_label_override
        if self.assembly_identifier:
            data['assembly'] = self.assembly_identifier
        return data

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemSelectedHasChanged:
            scene = self.scene()
            if isinstance(scene, SchematicScene):
                scene.on_item_selection_changed(self)
        return super().itemChange(change, value)

    def apply_data_label_visibility(self):
        visible = self.data_labels_visible_now()
        revealed_by_hover = self._hover_reveal and not self.data_labels_shown()
        for label in self.labels:
            label.setVisible(visible and bool(label.text))
            label.setZValue(20 if revealed_by_hover else 10)

    def set_points(self, points):
        self.points = [(float(point[0]), float(point[1])) for point in points]
        self.update_path()

    def update_path(self):
        self.prepareGeometryChange()
        path = QPainterPath()
        if self.points:
            path.moveTo(to_qpoint(self.points[0]))
            for point in self.points[1:]:
                path.lineTo(to_qpoint(point))
        self.setPath(path)
        self.layout_labels()
        self.setToolTip(f"<b>{self.name}</b> — inner diameter "
                        f"{format_number(float(self.properties.get('inner_diameter', 0)))} in, "
                        f"length {format_number(float(self.properties.get('length', 0)))} in, "
                        f"{int(self.properties.get('bend_count', 0))} bends<br>"
                        "Drag a segment to move it · drag an end to reconnect · right-click for menu")

    def shape(self):
        stroker = QPainterPathStroker()
        stroker.setWidth(12)
        return stroker.createStroke(self.path())

    def boundingRect(self):
        return self.path().boundingRect().adjusted(-16, -16, 16, 16)

    def set_result(self, result):
        self.result = result
        self.layout_labels()
        self.update()

    def remove_labels(self):
        for label in self.labels:
            if label.scene() is not None:
                label.scene().removeItem(label)
        self.labels = []

    def layout_labels(self):
        scene = self.scene()
        sections = (self.result or {}).get('sections', []) if len(self.points) >= 2 else []
        if scene is None:
            return
        while len(self.labels) < len(sections):
            label = LabelItem(None, point_size=7.5)
            scene.addItem(label)
            self.labels.append(label)
        while len(self.labels) > len(sections):
            label = self.labels.pop()
            if label.scene() is not None:
                label.scene().removeItem(label)
        self.arrows = []
        obstacles = []
        for block in getattr(scene, 'blocks', {}).values():
            obstacles.extend(block.label_rectangles())
        for other_line in getattr(scene, 'lines', {}).values():
            if other_line is self:
                continue
            for other_label in other_line.labels:
                if other_label.isVisible() and other_label.text:
                    center, rectangle = other_label.pos(), other_label.boundingRect()
                    obstacles.append((center.x() + rectangle.left(), center.y() + rectangle.top(),
                                      center.x() + rectangle.right(), center.y() + rectangle.bottom()))
        for label, section in zip(self.labels, sections):
            section_points = sub_polyline(self.points, section['start_fraction'], section['end_fraction'])
            longest_index = max(range(len(section_points) - 1),
                                key=lambda index: (abs(section_points[index + 1][0] - section_points[index][0])
                                                   + abs(section_points[index + 1][1] - section_points[index][1])))
            segment_start, segment_end = section_points[longest_index], section_points[longest_index + 1]
            middle = ((segment_start[0] + segment_end[0]) / 2, (segment_start[1] + segment_end[1]) / 2)
            mass_flow = section['mass_flow']
            if math.isnan(section['pressure_drop']):
                text = '— (no pressure reference)'
            elif abs(mass_flow) < 1e-9:
                text = '0 g/s'
            else:
                # pressure drop from inlet to outlet in the actual flow direction
                drop_in_flow_direction = (section['pressure_drop'] * (1 if mass_flow >= 0 else -1)
                                          / PASCALS_PER_PSI)
                text = (f"{format_number(abs(mass_flow) / KILOGRAMS_PER_GRAM, 3)} g/s\n"
                        f"ΔP {format_number(drop_in_flow_direction, 3)} psid")
            label.set_text(text)
            rectangle = label.boundingRect()
            center_x, center_y = place_label(segment_start, segment_end, rectangle.width(), rectangle.height(),
                                             obstacles)
            label.setPos(center_x, center_y)
            obstacles.append((center_x - rectangle.width() / 2, center_y - rectangle.height() / 2,
                              center_x + rectangle.width() / 2, center_y + rectangle.height() / 2))
            segment_length = abs(segment_end[0] - segment_start[0]) + abs(segment_end[1] - segment_start[1])
            if abs(mass_flow) > 1e-9 and segment_length > 0:
                flow_sign = 1 if mass_flow > 0 else -1
                self.arrows.append((middle, ((segment_end[0] - segment_start[0]) / segment_length * flow_sign,
                                             (segment_end[1] - segment_start[1]) / segment_length * flow_sign)))
        self.apply_data_label_visibility()
        self.update()

    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        body_color = (COLOR_LINE_SELECTED if self.isSelected()
                      else (QColor(THEME['gold_highlight']) if self._hovered else COLOR_LINE))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        # brass pipe: dark outline, gold body, light highlight
        for width, color in ((7.5, QColor(20, 12, 2, 235)), (4.2, body_color), (1.2, QColor(255, 244, 210, 150))):
            pen = QPen(color, width)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            painter.drawPath(self.path())
        if not self.points:
            return
        for attachment, end_point in ((self.start, self.points[0]), (self.end, self.points[-1])):
            if attachment and attachment.get('kind') == 'line':
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(body_color)
                painter.drawEllipse(to_qpoint(end_point), 5, 5)
            elif not attachment:
                painter.setPen(QPen(body_color, 2))
                painter.setBrush(QColor(THEME['background_deep']))
                painter.drawEllipse(to_qpoint(end_point), 4.5, 4.5)
        painter.setPen(QPen(QColor(4, 6, 12), 1.0))
        painter.setBrush(QColor(THEME['cyan']))
        for (middle_x, middle_y), (direction_x, direction_y) in self.arrows:
            normal_x, normal_y = -direction_y, direction_x
            painter.drawPolygon(QPolygonF([
                QPointF(middle_x + direction_x * 6, middle_y + direction_y * 6),
                QPointF(middle_x - direction_x * 5 + normal_x * 5, middle_y - direction_y * 5 + normal_y * 5),
                QPointF(middle_x - direction_x * 5 - normal_x * 5, middle_y - direction_y * 5 - normal_y * 5)]))
        if self.isSelected():
            painter.setPen(QPen(body_color, 1.2))
            painter.setBrush(QColor(THEME['background_deep']))
            for point in self.points[1:-1]:
                painter.drawRect(QRectF(point[0] - 3, point[1] - 3, 6, 6))

    def hoverEnterEvent(self, event):
        self._hovered = True
        self.set_hover_reveal(True)
        self.update()

    def hoverLeaveEvent(self, event):
        self._hovered = False
        self.set_hover_reveal(False)
        self.update()

    def hoverMoveEvent(self, event):
        position = to_tuple(event.pos())
        if self.points and (distance_between(position, self.points[0]) < 7
                            or distance_between(position, self.points[-1]) < 7):
            self.setCursor(Qt.CursorShape.CrossCursor)
        else:
            _, _, segment_index, _ = project_on_polyline(self.points, position)
            orientation = segment_orientation(self.points, segment_index) if len(self.points) > 1 else 'H'
            self.setCursor(Qt.CursorShape.SizeVerCursor if orientation == 'H' else Qt.CursorShape.SizeHorCursor)

    def mousePressEvent(self, event):
        scene = self.scene()
        self._drag = None
        if (event.button() == Qt.MouseButton.LeftButton and isinstance(scene, SchematicScene)
                and scene.mode == 'select' and len(self.points) >= 2 and not scene.is_locked_item(self)):
            position = to_tuple(event.pos())
            original_points = list(self.points)
            if distance_between(position, self.points[0]) < 7:
                self._drag = {'kind': 'end', 'which_end': 'start', 'original_points': original_points,
                              'moved': False}
            elif distance_between(position, self.points[-1]) < 7:
                self._drag = {'kind': 'end', 'which_end': 'end', 'original_points': original_points,
                              'moved': False}
            else:
                _, _, segment_index, _ = project_on_polyline(self.points, position)
                self._drag = {'kind': 'segment', 'segment_index': segment_index,
                              'original_points': original_points, 'press_position': position,
                              'orientation': segment_orientation(self.points, segment_index), 'moved': False}
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        scene = self.scene()
        drag = self._drag
        if not drag or not isinstance(scene, SchematicScene):
            return
        position = to_tuple(event.pos())
        if drag['kind'] == 'segment':
            offset = (snap_to_grid(position[1] - drag['press_position'][1]) if drag['orientation'] == 'H'
                      else snap_to_grid(position[0] - drag['press_position'][0]))
            new_points, _ = move_segment(drag['original_points'], drag['segment_index'], offset,
                                         start_fixed=True, end_fixed=True)
            if offset == 0:
                new_points = list(drag['original_points'])
            if not scene.collides(new_points):
                self.set_points(new_points)
                scene.update_dependents(self)
                drag['moved'] = drag['moved'] or offset != 0
        else:
            attachment, target_point, _ = scene.find_target(position, exclude_line=self)
            self.set_points(adjust_endpoint(drag['original_points'], drag['which_end'], target_point))
            drag['attachment'] = attachment
            drag['moved'] = True

    def mouseReleaseEvent(self, event):
        scene = self.scene()
        drag = self._drag
        self._drag = None
        super().mouseReleaseEvent(event)
        if not drag or not drag['moved'] or not isinstance(scene, SchematicScene):
            return
        if drag['kind'] == 'segment':
            self.set_points(tidy(self.points))
            self.auto_routed = False
            scene.update_dependents(self)
            scene.notify(physics_changed=False)
        else:
            attachment = drag.get('attachment')
            if attachment and attachment.get('kind') == 'line' and attachment.get('line') == self.identifier:
                attachment = None
            setattr(self, drag['which_end'], attachment)
            if self.auto_routed or scene.collides(self.points):
                scene.reroute_line(self)
            else:
                self.set_points(simplify(self.points))
            scene.refresh_port_markers()
            scene.notify()


# -----------------------------------------------------------------------------
# Scene
# -----------------------------------------------------------------------------
class SchematicScene(QGraphicsScene):
    model_changed = Signal(bool)    # True when the physics may have changed
    message = Signal(str)

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.setSceneRect(-3000, -2000, 9000, 6000)
        self.blocks = {}
        self.lines = {}
        self.mode = 'select'
        self.connection = None      # port drag in progress
        self.drawing = None         # draw-mode polyline in progress
        self.show_labels = True
        self.fluid = dict(DEFAULT_FLUID)
        self.settings = dict(DEFAULT_SETTINGS)
        self.scenarios = []
        self.active_scenario = None
        self.loading = False
        self.results = None
        self._moved_block_identifiers = set()
        self.assemblies = {}            # identifier -> {'name': str, 'locked': bool}
        self.assembly_outlines = {}     # identifier -> AssemblyOutlineItem
        self._expanding_selection = False
        self.preview = QGraphicsPathItem()
        self.preview.setZValue(50)
        self.addItem(self.preview)
        self.preview.hide()

    # ---- infrastructure ----------------------------------------------------
    def notify(self, physics_changed=True):
        if not self.loading:
            self.update_assembly_outlines()
            self.model_changed.emit(physics_changed)

    # ---- assemblies: saved groups of connected parts, locked or unlocked -------------
    def assembly_record(self, item):
        return self.assemblies.get(getattr(item, 'assembly_identifier', None))

    def is_locked_item(self, item):
        record = self.assembly_record(item)
        return bool(record and record.get('locked'))

    def assembly_members(self, assembly_identifier):
        return ([block for block in self.blocks.values() if block.assembly_identifier == assembly_identifier],
                [line for line in self.lines.values() if line.assembly_identifier == assembly_identifier])

    def on_item_selection_changed(self, item):
        """Selecting any part of a locked assembly selects the whole assembly."""
        if self._expanding_selection or self.loading or not item.isSelected() or not self.is_locked_item(item):
            return
        self._expanding_selection = True
        try:
            member_blocks, member_lines = self.assembly_members(item.assembly_identifier)
            for member in member_blocks + member_lines:
                if not member.isSelected():
                    member.setSelected(True)
        finally:
            self._expanding_selection = False

    def create_assembly(self, items, name, locked=True):
        items = [item for item in items if isinstance(item, (BlockItem, LineItem))]
        if not items:
            return None
        assembly_identifier = new_identifier()
        self.assemblies[assembly_identifier] = {'name': name, 'locked': bool(locked)}
        for item in items:
            item.assembly_identifier = assembly_identifier
        self._prune_assemblies()
        self.notify(physics_changed=False)
        return assembly_identifier

    def set_assembly_locked(self, assembly_identifier, locked):
        if assembly_identifier in self.assemblies:
            self.assemblies[assembly_identifier]['locked'] = bool(locked)
            self.notify(physics_changed=False)

    def rename_assembly(self, assembly_identifier, name):
        if assembly_identifier in self.assemblies and name:
            self.assemblies[assembly_identifier]['name'] = name
            self.notify(physics_changed=False)

    def add_items_to_assembly(self, assembly_identifier, items):
        for item in items:
            if isinstance(item, (BlockItem, LineItem)):
                item.assembly_identifier = assembly_identifier
        self._prune_assemblies()
        self.notify(physics_changed=False)

    def remove_items_from_assembly(self, items):
        for item in items:
            if isinstance(item, (BlockItem, LineItem)):
                item.assembly_identifier = None
        self._prune_assemblies()
        self.notify(physics_changed=False)

    def dissolve_assembly(self, assembly_identifier):
        member_blocks, member_lines = self.assembly_members(assembly_identifier)
        self.remove_items_from_assembly(member_blocks + member_lines)

    def _prune_assemblies(self):
        """Forget assemblies that no longer have any parts."""
        used = {item.assembly_identifier for item in list(self.blocks.values()) + list(self.lines.values())}
        for assembly_identifier in list(self.assemblies):
            if assembly_identifier not in used:
                del self.assemblies[assembly_identifier]

    def update_assembly_outlines(self):
        for assembly_identifier in list(self.assembly_outlines):
            if assembly_identifier not in self.assemblies:
                self.removeItem(self.assembly_outlines.pop(assembly_identifier))
        for assembly_identifier, record in self.assemblies.items():
            member_blocks, member_lines = self.assembly_members(assembly_identifier)
            rectangles = [rectangle for block in member_blocks for rectangle in block.label_rectangles()]
            for line in member_lines:
                if line.points:
                    xs = [point[0] for point in line.points]
                    ys = [point[1] for point in line.points]
                    rectangles.append((min(xs), min(ys), max(xs), max(ys)))
            if not rectangles:
                continue
            padding = 12.0
            left = min(rectangle[0] for rectangle in rectangles) - padding
            top = min(rectangle[1] for rectangle in rectangles) - padding
            right = max(rectangle[2] for rectangle in rectangles) + padding
            bottom = max(rectangle[3] for rectangle in rectangles) + padding
            outline_rectangle = QRectF(left, top, right - left, bottom - top)
            outline = self.assembly_outlines.get(assembly_identifier)
            if outline is None:
                outline = AssemblyOutlineItem(assembly_identifier)
                self.addItem(outline)
                self.assembly_outlines[assembly_identifier] = outline
            outline.set_geometry(outline_rectangle, record['name'], record['locked'])

    def drawBackground(self, painter, rectangle):
        background = QLinearGradient(QPointF(0, rectangle.top()), QPointF(0, rectangle.bottom()))
        background.setColorAt(0.0, QColor('#060a16'))
        background.setColorAt(1.0, QColor(THEME['background_deep']))
        painter.fillRect(rectangle, QBrush(background))
        grid_step = 20
        if (rectangle.width() / grid_step) * (rectangle.height() / grid_step) > 60000:
            return
        painter.setPen(QPen(QColor(70, 150, 190, 70), 1.6))
        first_x = int(math.floor(rectangle.left() / grid_step)) * grid_step
        first_y = int(math.floor(rectangle.top() / grid_step)) * grid_step
        grid_points = [QPointF(x, y) for x in range(first_x, int(rectangle.right()) + 1, grid_step)
                       for y in range(first_y, int(rectangle.bottom()) + 1, grid_step)]
        painter.drawPoints(QPolygonF(grid_points))
        painter.setPen(QPen(QColor(217, 164, 65, 60), 1.0))
        crosshair_step = 100
        for x in range(int(math.floor(rectangle.left() / crosshair_step)) * crosshair_step,
                       int(rectangle.right()) + 1, crosshair_step):
            for y in range(int(math.floor(rectangle.top() / crosshair_step)) * crosshair_step,
                           int(rectangle.bottom()) + 1, crosshair_step):
                painter.drawLine(QPointF(x - 4, y), QPointF(x + 4, y))
                painter.drawLine(QPointF(x, y - 4), QPointF(x, y + 4))

    def _set_preview(self, points, valid=True):
        path = QPainterPath()
        if points:
            path.moveTo(to_qpoint(points[0]))
            for point in points[1:]:
                path.lineTo(to_qpoint(point))
        self.preview.setPath(path)
        self.preview.setPen(QPen(QColor(THEME['cyan']) if valid else QColor(THEME['warning']), 2.5,
                                 Qt.PenStyle.DashLine))
        self.preview.show()

    def set_mode(self, mode):
        self.cancel_drawing()
        self.connection = None
        self.preview.hide()
        self.mode = mode
        for view in self.views():
            view.setDragMode(QGraphicsView.DragMode.RubberBandDrag if mode == 'select'
                             else QGraphicsView.DragMode.NoDrag)
            view.viewport().setCursor(Qt.CursorShape.ArrowCursor if mode == 'select'
                                      else Qt.CursorShape.CrossCursor)
        self.message.emit("Draw mode: click to place horizontal/vertical sections, click a port or "
                          "line to finish, right-click/double-click to end free, Esc to cancel."
                          if mode == 'draw' else "")

    # ---- geometry helpers ----------------------------------------------------
    def block_rectangles(self):
        return [block.scene_rectangle() for block in self.blocks.values()]

    def port_positions(self):
        return [block.port_scene_position(port_name) for block in self.blocks.values() for port_name in block.ports]

    def collides(self, points):
        """True if the polyline passes over a block body or a port it is not attached to."""
        if len(points) < 2:
            return False
        foreign_ports = [position for position in self.port_positions()
                         if not points_close(position, points[0]) and not points_close(position, points[-1])]
        return path_crosses_obstacles(points, self.block_rectangles(), foreign_ports)

    def route(self, start, start_direction, end, end_direction):
        foreign_ports = [position for position in self.port_positions()
                         if not points_close(position, start) and not points_close(position, end)]
        points = route_around_blocks(start, start_direction, end, end_direction, self.block_rectangles(),
                                     foreign_ports)
        if points is None:
            self.message.emit("No clear route found — line drawn directly; drag its segments clear.")
            points = quick_orthogonal_route(start, start_direction, end, end_direction)
        return points

    def attachment_point(self, attachment):
        """Scene position and outward direction of an attachment (port or line tap)."""
        if not attachment:
            return None, None
        if attachment.get('kind') == 'port':
            block = self.blocks.get(attachment.get('block'))
            if block and attachment.get('port') in block.ports:
                return block.port_scene_position(attachment['port']), block.port_direction(attachment['port'])
        if attachment.get('kind') == 'line':
            line = self.lines.get(attachment.get('line'))
            if line and len(line.points) >= 2:
                return point_at_fraction(line.points, float(attachment.get('path_fraction', 0)))[0], None
        return None, None

    def reroute_line(self, line):
        start, start_direction = self.attachment_point(line.start)
        end, end_direction = self.attachment_point(line.end)
        start = start or line.points[0]
        end = end or line.points[-1]
        line.set_points(self.route(start, start_direction, end, end_direction))
        self.update_dependents(line)

    def update_dependents(self, line, visited=None):
        """Move the ends of lines tapped onto `line` so they follow it."""
        visited = visited if visited is not None else set()
        visited.add(line.identifier)
        for other_line in list(self.lines.values()):
            if other_line.identifier in visited or len(other_line.points) < 2:
                continue
            changed = False
            for which_end in ('start', 'end'):
                attachment = getattr(other_line, which_end)
                if attachment and attachment.get('kind') == 'line' and attachment.get('line') == line.identifier:
                    tap_point, _ = point_at_fraction(line.points, float(attachment.get('path_fraction', 0)))
                    current = other_line.points[0] if which_end == 'start' else other_line.points[-1]
                    if not points_close(tap_point, current):
                        other_line.points = adjust_endpoint(other_line.points, which_end, tap_point)
                        changed = True
            if changed:
                other_line.update_path()
                self.update_dependents(other_line, visited)

    def find_target(self, position, exclude_line=None, exclude_port=None):
        """What is under `position` for connecting: (attachment or None, point, direction)."""
        best = None
        for block in self.blocks.values():
            for port_name, port in block.ports.items():
                if port is exclude_port:
                    continue
                port_position = block.port_scene_position(port_name)
                distance = distance_between(port_position, position)
                if distance <= 10 and (best is None or distance < best[0]):
                    best = (distance, {'kind': 'port', 'block': block.identifier, 'port': port_name},
                            port_position, block.port_direction(port_name))
        if best:
            return best[1], best[2], best[3]
        for line in self.lines.values():
            if line is exclude_line or len(line.points) < 2:
                continue
            path_fraction, closest, segment_index, distance = project_on_polyline(line.points, position)
            if distance <= 6:
                segment_start, segment_end = line.points[segment_index], line.points[segment_index + 1]
                if segment_orientation(line.points, segment_index) == 'H':
                    closest = (min(max(snap_to_grid(closest[0]), min(segment_start[0], segment_end[0])),
                                   max(segment_start[0], segment_end[0])), segment_start[1])
                else:
                    closest = (segment_start[0], min(max(snap_to_grid(closest[1]), min(segment_start[1], segment_end[1])),
                                                     max(segment_start[1], segment_end[1])))
                path_fraction = project_on_polyline(line.points, closest)[0]
                path_fraction = 0.0 if path_fraction < 1e-6 else (1.0 if path_fraction > 1 - 1e-6 else path_fraction)
                return {'kind': 'line', 'line': line.identifier, 'path_fraction': path_fraction}, closest, None
        return None, (snap_to_grid(position[0]), snap_to_grid(position[1])), None

    def tubing_properties_for(self, attachment):
        """Outer / inner diameter for a new line, taken from the port or line it starts on."""
        if attachment and attachment.get('kind') == 'port':
            block = self.blocks.get(attachment.get('block'))
            tubing_key = PORT_TUBING_PROPERTY.get(block.block_type, {}).get(attachment.get('port')) if block else None
            if tubing_key:
                outer_diameter, wall_thickness = tube_outer_diameter_and_wall(
                    block.properties.get(tubing_key, DEFAULT_TUBE))
                return {'outer_diameter': outer_diameter, 'inner_diameter': outer_diameter - 2 * wall_thickness}
        if attachment and attachment.get('kind') == 'line':
            line = self.lines.get(attachment.get('line'))
            if line:
                return {'outer_diameter': line.properties.get('outer_diameter'),
                        'inner_diameter': line.properties.get('inner_diameter')}
        return {}

    def refresh_port_markers(self):
        connected_ports = set()
        for line in self.lines.values():
            for attachment in (line.start, line.end):
                if attachment and attachment.get('kind') == 'port':
                    connected_ports.add((attachment.get('block'), attachment.get('port')))
        for block in self.blocks.values():
            for port_name, port in block.ports.items():
                port.set_connected((block.identifier, port_name) in connected_ports)

    # ---- block movement -----------------------------------------------------
    def on_block_moving(self, block, delta=(0.0, 0.0)):
        if self.loading:
            return
        self._moved_block_identifiers.add(block.identifier)
        locked_assembly = block.assembly_identifier if self.is_locked_item(block) else None
        if locked_assembly and (delta[0] or delta[1]):
            member_blocks, member_lines = self.assembly_members(locked_assembly)
            # a locked assembly moves rigidly: its own tubes are carried along unchanged
            # (done once per step, driven by the member with the lowest identifier)
            if block.identifier == min(member.identifier for member in member_blocks):
                for line in member_lines:
                    line.points = [(point[0] + delta[0], point[1] + delta[1]) for point in line.points]
                    line.update_path()
                    self.update_dependents(line)
        for line in self.lines.values():
            if locked_assembly and line.assembly_identifier == locked_assembly:
                continue
            changed = False
            for which_end in ('start', 'end'):
                attachment = getattr(line, which_end)
                if (attachment and attachment.get('kind') == 'port' and attachment.get('block') == block.identifier
                        and len(line.points) >= 2):
                    line.points = adjust_endpoint(line.points, which_end,
                                                  block.port_scene_position(attachment['port']))
                    changed = True
            if changed:
                line.update_path()
                self.update_dependents(line)

    def on_block_move_finished(self):
        moved, self._moved_block_identifiers = self._moved_block_identifiers, set()
        if not moved or self.loading:
            return
        for line in list(self.lines.values()):
            if len(line.points) < 2:
                continue
            attached = any(attachment and attachment.get('kind') == 'port' and attachment.get('block') in moved
                           for attachment in (line.start, line.end))
            carried = self.is_locked_item(line)       # moved rigidly with its locked assembly
            if self.collides(line.points) or (attached and line.auto_routed and not carried):
                self.reroute_line(line)
            elif attached and not carried:
                line.set_points(simplify(line.points))
        self.notify(physics_changed=False)

    # ---- creation / deletion ------------------------------------------------
    def next_name(self, prefix, for_lines=False):
        names = ([line.name for line in self.lines.values()] if for_lines
                 else [block.name for block in self.blocks.values()])
        pattern = re.compile(rf'^{re.escape(prefix)}-?(\d+)$')
        highest = 0
        for name in names:
            match = pattern.match(name or '')
            if match:
                highest = max(highest, int(match.group(1)))
        return f"{prefix}{highest + 1}" if for_lines else f"{prefix}-{highest + 1}"

    def add_block(self, data, position=None, select=True):
        block_data = copy.deepcopy(data)
        if not block_data.get('identifier') or block_data['identifier'] in self.blocks:
            block_data['identifier'] = new_identifier()
        if not block_data.get('name'):
            block_data['name'] = self.next_name(BLOCK_TYPES[block_data['type']]['prefix'])
        if position is not None:
            block_data['x'], block_data['y'] = snap_to_grid(position[0]), snap_to_grid(position[1])
        block = BlockItem(block_data)
        self.addItem(block)
        self.blocks[block.identifier] = block
        block.set_result(None)
        if not self.loading:
            if select:
                self.clearSelection()
                block.setSelected(True)
            self._moved_block_identifiers.add(block.identifier)
            self.on_block_move_finished()
            self.refresh_port_markers()
            self.notify()
        return block

    def add_line(self, data):
        line_data = copy.deepcopy(data)
        if not line_data.get('identifier') or line_data['identifier'] in self.lines:
            line_data['identifier'] = new_identifier()
        if not line_data.get('name'):
            line_data['name'] = self.next_name('L', for_lines=True)
        line = LineItem(line_data)
        self.addItem(line)
        self.lines[line.identifier] = line
        line.update_path()
        if not self.loading:
            self.refresh_port_markers()
        return line

    def remove_block(self, block):
        for line in self.lines.values():
            for which_end in ('start', 'end'):
                attachment = getattr(line, which_end)
                if attachment and attachment.get('kind') == 'port' and attachment.get('block') == block.identifier:
                    setattr(line, which_end, None)
                    line.update()
        self.removeItem(block)
        self.blocks.pop(block.identifier, None)

    def remove_line(self, line):
        for other_line in self.lines.values():
            for which_end in ('start', 'end'):
                attachment = getattr(other_line, which_end)
                if attachment and attachment.get('kind') == 'line' and attachment.get('line') == line.identifier:
                    setattr(other_line, which_end, None)
                    other_line.update()
        line.remove_labels()
        self.removeItem(line)
        self.lines.pop(line.identifier, None)

    def delete_items(self, items):
        for item in list(items):
            if isinstance(item, BlockItem) and item.identifier in self.blocks:
                self.remove_block(item)
            elif isinstance(item, LineItem) and item.identifier in self.lines:
                self.remove_line(item)
        self._prune_assemblies()
        self.refresh_port_markers()
        self.notify()

    # ---- copy / paste -------------------------------------------------------
    def copy_payload(self, items):
        """Clipboard dictionary for the given blocks and lines.

        Lines come along when they are selected, or when every attached end lands on a
        copied block or another copied line (so connecting tubes and transducer taps are
        copied with a group).  An end attached to something outside the copy becomes free."""
        blocks = [item for item in items if isinstance(item, BlockItem)]
        block_identifiers = {block.identifier for block in blocks}
        line_identifiers = {item.identifier for item in items if isinstance(item, LineItem)}

        def end_inside(attachment):
            if not attachment:
                return None                                   # free end: neutral
            if attachment.get('kind') == 'port':
                return attachment.get('block') in block_identifiers
            if attachment.get('kind') == 'line':
                return attachment.get('line') in line_identifiers
            return False

        added = True
        while added:                                          # grow until no more lines qualify
            added = False
            for line in self.lines.values():
                if line.identifier in line_identifiers:
                    continue
                ends = [end_inside(line.start), end_inside(line.end)]
                touches_block = any(attachment and attachment.get('kind') == 'port'
                                    and attachment.get('block') in block_identifiers
                                    for attachment in (line.start, line.end))
                if touches_block and all(inside is not False for inside in ends):
                    line_identifiers.add(line.identifier)
                    added = True
        line_dictionaries = []
        for line in self.lines.values():
            if line.identifier not in line_identifiers:
                continue
            line_dictionary = line.to_dictionary()
            for which_end in ('start', 'end'):
                if end_inside(line_dictionary[which_end]) is False:
                    line_dictionary[which_end] = None
                    line_dictionary['auto_routed'] = False
            line_dictionaries.append(line_dictionary)
        if not blocks and not line_dictionaries:
            return None
        block_dictionaries = [block.to_dictionary() for block in blocks]
        referenced = {item.get('assembly') for item in block_dictionaries + line_dictionaries} - {None}
        return {'format': 'pylon_clipboard', 'version': FILE_FORMAT_VERSION,
                'blocks': block_dictionaries, 'lines': line_dictionaries,
                'assemblies': {identifier: dict(self.assemblies[identifier])
                               for identifier in referenced if identifier in self.assemblies}}

    @staticmethod
    def _payload_bounds(payload):
        """(left, top, right, bottom) of everything in a clipboard payload."""
        rectangles = [block_rectangle(block) for block in payload['blocks']]
        for line in payload['lines']:
            xs = [point[0] for point in line['points']]
            ys = [point[1] for point in line['points']]
            if xs:
                rectangles.append((min(xs), min(ys), max(xs), max(ys)))
        return (min(rectangle[0] for rectangle in rectangles), min(rectangle[1] for rectangle in rectangles),
                max(rectangle[2] for rectangle in rectangles), max(rectangle[3] for rectangle in rectangles))

    def _offset_is_clear(self, payload, offset_x, offset_y):
        """True if the payload moved by the offset neither overlaps a block nor crosses a line/block."""
        clearance = GRID * 2
        existing_rectangles = self.block_rectangles()
        moved_rectangles = []
        for block in payload['blocks']:
            left, top, right, bottom = block_rectangle(block)
            moved_rectangles.append((left + offset_x, top + offset_y, right + offset_x, bottom + offset_y))
        for left, top, right, bottom in moved_rectangles:
            for other_left, other_top, other_right, other_bottom in existing_rectangles:
                if (left - clearance < other_right and other_left < right + clearance
                        and top - clearance < other_bottom and other_top < bottom + clearance):
                    return False
        existing_ports = self.port_positions()
        for line in payload['lines']:
            moved_points = [(point[0] + offset_x, point[1] + offset_y) for point in line['points']]
            if len(moved_points) >= 2 and path_crosses_obstacles(moved_points, existing_rectangles, existing_ports):
                return False
        for line in self.lines.values():
            if len(line.points) >= 2 and path_crosses_obstacles(line.points, moved_rectangles):
                return False
        return True

    def paste_payload(self, payload, position=None, step=1):
        """Add a copy of a clipboard payload; returns the new items.

        With `position` the group is centred there; otherwise it is offset diagonally by
        `step` grid steps of 40 px from where it was copied.  Either way it is then nudged
        right/down until it sits clear of existing blocks and lines."""
        if not payload or not (payload.get('blocks') or payload.get('lines')):
            return []
        source_assemblies = dict(payload.get('assemblies') or {})
        payload = {'blocks': [upgrade_block_template(block) for block in payload.get('blocks', [])
                              if block.get('type') in BLOCK_TYPES],
                   'lines': copy.deepcopy(payload.get('lines', []))}
        if not payload['blocks'] and not payload['lines']:
            return []
        left, top, right, bottom = self._payload_bounds(payload)
        if position is not None:
            base_x = snap_to_grid(position[0] - (left + right) / 2)
            base_y = snap_to_grid(position[1] - (top + bottom) / 2)
        else:
            base_x = base_y = 40 * step
        offset_x, offset_y = base_x, base_y
        for attempt in range(1, 80):
            if self._offset_is_clear(payload, offset_x, offset_y):
                break
            ring, side = divmod(attempt, 3)
            nudge = 40 * (ring + 1)
            offset_x, offset_y = [(base_x + nudge, base_y), (base_x, base_y + nudge),
                                  (base_x + nudge, base_y + nudge)][side]
        else:
            # no clear spot nearby: place the group below everything already drawn
            lowest = max([rectangle[3] for rectangle in self.block_rectangles()]
                         + [point[1] for line in self.lines.values() for point in line.points] + [bottom])
            offset_x, offset_y = base_x, snap_to_grid(lowest + 80 - top)
        block_identifier_map = {block['identifier']: new_identifier() for block in payload['blocks']}
        line_identifier_map = {line['identifier']: new_identifier() for line in payload['lines']}

        def remap(attachment):
            if not attachment:
                return None
            attachment = dict(attachment)
            if attachment.get('kind') == 'port':
                if attachment.get('block') not in block_identifier_map:
                    return None
                attachment['block'] = block_identifier_map[attachment['block']]
            elif attachment.get('kind') == 'line':
                if attachment.get('line') not in line_identifier_map:
                    return None
                attachment['line'] = line_identifier_map[attachment['line']]
            return attachment

        # pasted assembly parts form new assemblies (names incremented like block names)
        taken_assembly_names = {record['name'] for record in self.assemblies.values()}
        assembly_identifier_map = {}
        for old_identifier, record in source_assemblies.items():
            assembly_name = record.get('name') or 'Assembly'
            if assembly_name in taken_assembly_names:
                assembly_name = incremented_name(assembly_name, taken_assembly_names)
            taken_assembly_names.add(assembly_name)
            assembly_identifier_map[old_identifier] = new_identifier()
            self.assemblies[assembly_identifier_map[old_identifier]] = {'name': assembly_name,
                                                                        'locked': bool(record.get('locked', True))}
        new_items = []
        taken_names = {existing.name for existing in self.blocks.values()}
        self.loading = True                       # defer re-routing until everything exists
        try:
            for block in payload['blocks']:
                block_data = copy.deepcopy(block)
                block_data['identifier'] = block_identifier_map[block['identifier']]
                # copy name: last number + 1 (or ' 1' added), skipping names already in use
                block_data['name'] = (incremented_name(block['name'], taken_names) if block.get('name')
                                      else None)
                if block_data['name']:
                    taken_names.add(block_data['name'])
                block_data['assembly'] = assembly_identifier_map.get(block.get('assembly'))
                block_data['x'] = block['x'] + offset_x
                block_data['y'] = block['y'] + offset_y
                new_items.append(self.add_block(block_data))
            for line in payload['lines']:
                line_data = copy.deepcopy(line)
                line_data['identifier'] = line_identifier_map[line['identifier']]
                line_data['assembly'] = assembly_identifier_map.get(line.get('assembly'))
                line_data['name'] = None
                line_data['points'] = [(point[0] + offset_x, point[1] + offset_y) for point in line['points']]
                line_data['start'] = remap(line['start'])
                line_data['end'] = remap(line['end'])
                new_items.append(self.add_line(line_data))
        finally:
            self.loading = False
        for line in list(self.lines.values()):     # nothing may pass over a block
            if len(line.points) >= 2 and self.collides(line.points):
                self.reroute_line(line)
        self.clearSelection()
        for item in new_items:
            item.setSelected(True)
        self.refresh_port_markers()
        self.notify()
        return new_items

    def duplicate_items(self, items):
        """Duplicate blocks (and their connecting lines) next to the originals."""
        payload = self.copy_payload(items)
        if payload:
            self.paste_payload(payload)

    def duplicate_blocks(self, blocks):
        self.duplicate_items(blocks)

    def rotate_block(self, block, clockwise=True):
        """Rotate 90° to the right (clockwise on screen) or to the left."""
        if self.is_locked_item(block):
            self.message.emit(f"'{self.assembly_record(block)['name']}' is locked — unlock the assembly to "
                              "rotate its parts.")
            return
        block.orientation = (block.orientation + (90 if clockwise else 270)) % 360
        block.layout_ports()
        block.update()
        self.on_block_moving(block)
        self.on_block_move_finished()

    # ---- port-drag connections ---------------------------------------------
    def begin_connection(self, port):
        self.connection = {'port': port}

    def update_connection(self, scene_position):
        port = self.connection['port']
        start = port.block.port_scene_position(port.name)
        attachment, target_point, target_direction = self.find_target(to_tuple(scene_position), exclude_port=port)
        self._set_preview(quick_orthogonal_route(start, port.direction, target_point, target_direction),
                          valid=attachment is not None)

    def end_connection(self, scene_position):
        self.preview.hide()
        connection, self.connection = self.connection, None
        if not connection:
            return
        port = connection['port']
        block = port.block
        start = block.port_scene_position(port.name)
        attachment, target_point, target_direction = self.find_target(to_tuple(scene_position), exclude_port=port)
        if attachment is None and distance_between(start, target_point) < 15:
            return
        start_attachment = {'kind': 'port', 'block': block.identifier, 'port': port.name}
        properties = default_line_properties()
        properties.update(self.tubing_properties_for(start_attachment))
        self.add_line({'points': self.route(start, port.direction, target_point, target_direction),
                       'start': start_attachment, 'end': attachment, 'properties': properties,
                       'auto_routed': True})
        self.notify()

    # ---- draw mode -----------------------------------------------------------
    @staticmethod
    def _orthogonal_projection(last_point, position):
        if abs(position[0] - last_point[0]) >= abs(position[1] - last_point[1]):
            return (snap_to_grid(position[0]), last_point[1])
        return (last_point[0], snap_to_grid(position[1]))

    @staticmethod
    def _previous_orientation(points):
        if len(points) >= 2:
            return segment_orientation(points, len(points) - 2)
        return None

    def drawing_click(self, position):
        attachment, target_point, _ = self.find_target(position)
        if self.drawing is None:
            self.drawing = {'points': [target_point], 'start': attachment}
            self._set_preview([target_point, target_point])
            return
        points = self.drawing['points']
        last_point = points[-1]
        if attachment is not None and not points_close(target_point, points[0]):
            new_points = points + orthogonal_connection(last_point, target_point, self._previous_orientation(points))
            if self.collides(new_points):
                self.message.emit("A line can't pass over a component — add a bend first.")
                return
            self._finish_drawing(new_points, attachment)
            return
        next_point = self._orthogonal_projection(last_point, position)
        if points_close(next_point, last_point):
            return
        new_points = points + [next_point]
        if self.collides(new_points):
            self.message.emit("A line can't pass over a component — choose another path.")
            return
        self.drawing['points'] = new_points

    def drawing_move(self, position):
        if not self.drawing:
            return
        points = self.drawing['points']
        attachment, target_point, _ = self.find_target(position)
        last_point = points[-1]
        tail = (orthogonal_connection(last_point, target_point, self._previous_orientation(points)) if attachment
                else [self._orthogonal_projection(last_point, position)])
        preview_points = points + tail
        self._set_preview(preview_points, valid=not self.collides(preview_points))

    def finish_drawing_with_free_end(self):
        if self.drawing and len(self.drawing['points']) >= 2:
            self._finish_drawing(self.drawing['points'], None)
        else:
            self.cancel_drawing()

    def _finish_drawing(self, points, end_attachment):
        start_attachment = self.drawing['start']
        self.drawing = None
        self.preview.hide()
        points = simplify(points)
        if len(points) < 2 or polyline_length(points) == 0:
            return
        properties = default_line_properties()
        properties.update(self.tubing_properties_for(start_attachment)
                          or self.tubing_properties_for(end_attachment))
        self.add_line({'points': points, 'start': start_attachment, 'end': end_attachment,
                       'properties': properties, 'auto_routed': False})
        self.notify()

    def cancel_drawing(self):
        self.drawing = None
        self.preview.hide()

    def escape(self):
        if self.drawing:
            self.cancel_drawing()
        elif self.mode == 'draw':
            self.main_window.set_mode('select')
        else:
            self.clearSelection()

    # ---- line editing ---------------------------------------------------------
    def add_segment(self, line, segment_index, position):
        line.set_points(split_segment(line.points, segment_index, position))
        line.auto_routed = False
        self.message.emit("Segment split — drag either half to create the new section.")
        self.notify(physics_changed=False)

    def delete_line_segment(self, line, segment_index):
        new_points = delete_segment(line.points, segment_index, line.start is not None, line.end is not None)
        if new_points is None:
            self.message.emit("That segment can't be removed without moving an attached end.")
            return
        if self.collides(new_points):
            self.message.emit("Removing that segment would route the line over a component.")
            return
        line.set_points(new_points)
        line.auto_routed = False
        self.update_dependents(line)
        self.notify(physics_changed=False)

    def auto_route(self, lines):
        for line in lines:
            line.auto_routed = True
            self.reroute_line(line)
        self.notify(physics_changed=False)

    def reverse_line(self, line):
        line.points = list(reversed(line.points))
        line.start, line.end = line.end, line.start
        for other_line in self.lines.values():
            for which_end in ('start', 'end'):
                attachment = getattr(other_line, which_end)
                if attachment and attachment.get('kind') == 'line' and attachment.get('line') == line.identifier:
                    attachment['path_fraction'] = 1.0 - float(attachment.get('path_fraction', 0))
        line.update_path()
        self.notify()

    def combine_lines(self, selected_lines):
        """Join lines that form one end-to-end chain into a single line."""
        line_ends = []
        for line in selected_lines:
            line_ends.append((line, 'start', line.points[0]))
            line_ends.append((line, 'end', line.points[-1]))

        def partner_of(line, point):
            for other_line, other_end, other_point in line_ends:
                if other_line is not line and points_close(other_point, point, 1.0):
                    return other_line, other_end
            return None
        free_ends = [(line, which_end) for line, which_end, point in line_ends if partner_of(line, point) is None]
        if len(free_ends) != 2:
            QMessageBox.information(self.main_window, "Combine lines",
                                    "The selected lines must form a single continuous chain "
                                    "(end-to-end, no branches).")
            return
        current_line, entry_end = free_ends[0]
        used_identifiers, chain_order, combined_points = set(), [], []
        while current_line is not None and current_line.identifier not in used_identifiers:
            used_identifiers.add(current_line.identifier)
            sequence = (list(current_line.points) if entry_end == 'start'
                        else list(reversed(current_line.points)))
            combined_points = combined_points + (sequence[1:] if combined_points else sequence)
            chain_order.append((current_line, entry_end))
            next_link = partner_of(current_line, sequence[-1])
            if next_link is None:
                break
            current_line, entry_end = next_link
        if len(used_identifiers) != len(selected_lines):
            QMessageBox.information(self.main_window, "Combine lines", "The selected lines do not form one chain.")
            return
        first_line, first_entry = chain_order[0]
        last_line, last_entry = chain_order[-1]
        start_attachment = first_line.start if first_entry == 'start' else first_line.end
        end_attachment = last_line.end if last_entry == 'start' else last_line.start
        properties = dict(first_line.properties)
        properties['length'] = sum(float(line.properties.get('length', 0)) for line in selected_lines)
        properties['bend_count'] = int(sum(float(line.properties.get('bend_count', 0)) for line in selected_lines))
        dialog = PropertyDialog(self.main_window, "Combined line properties", LINE_PROPERTIES, properties,
                                name=first_line.name, tube_helper=True,
                                note="Length and bends default to the sums of the selected lines.")
        if not dialog.exec():
            return
        if start_attachment and start_attachment.get('kind') == 'line' and start_attachment.get('line') in used_identifiers:
            start_attachment = None
        if end_attachment and end_attachment.get('kind') == 'line' and end_attachment.get('line') in used_identifiers:
            end_attachment = None
        combined_line = self.add_line({'points': simplify(combined_points), 'start': copy.deepcopy(start_attachment),
                                       'end': copy.deepcopy(end_attachment), 'properties': dialog.values,
                                       'auto_routed': False, 'name': dialog.name})
        for other_line in self.lines.values():
            if other_line is combined_line or other_line.identifier in used_identifiers:
                continue
            for which_end in ('start', 'end'):
                attachment = getattr(other_line, which_end)
                if attachment and attachment.get('kind') == 'line' and attachment.get('line') in used_identifiers:
                    end_point = other_line.points[0] if which_end == 'start' else other_line.points[-1]
                    setattr(other_line, which_end, {'kind': 'line', 'line': combined_line.identifier,
                                                    'path_fraction': project_on_polyline(combined_line.points,
                                                                                         end_point)[0]})
        for line in selected_lines:
            line.remove_labels()
            self.removeItem(line)
            self.lines.pop(line.identifier, None)
        self.clearSelection()
        combined_line.setSelected(True)
        self.refresh_port_markers()
        self.notify()

    # ---- context menus --------------------------------------------------------
    def _target_at(self, scene_position):
        for item in self.items(scene_position):
            current = item
            while current is not None:
                if isinstance(current, (BlockItem, LineItem)):
                    return current
                current = current.parentItem()
        return None

    def contextMenuEvent(self, event):
        if self.mode == 'draw':
            event.accept()
            return
        target = self._target_at(event.scenePos())
        if isinstance(target, BlockItem):
            self.block_menu(target, event)
        elif isinstance(target, LineItem):
            self.line_menu(target, event)
        else:
            self.canvas_menu(event)
        event.accept()

    def _add_data_label_action(self, menu, clicked_item):
        """Checkable 'Show flow / pressure data' for the clicked element (and the rest of the selection)."""
        targets = [item for item in self.selectedItems() if isinstance(item, (BlockItem, LineItem))]
        if clicked_item not in targets:
            targets = [clicked_item]
        action = menu.addAction("Show flow / pressure data here")
        action.setCheckable(True)
        action.setChecked(clicked_item.data_labels_shown())

        def apply(visible, targets=targets):
            for target in targets:
                target.set_data_label_override(visible)
            self.notify(physics_changed=False)
        action.triggered.connect(apply)

    def _add_assembly_actions(self, menu, clicked_item):
        """Assembly submenu for a part of an assembly, or 'Save selection as assembly…'."""
        selected = [item for item in self.selectedItems() if isinstance(item, (BlockItem, LineItem))] or [clicked_item]
        assembly_identifier = clicked_item.assembly_identifier
        record = self.assemblies.get(assembly_identifier)
        menu.addSeparator()
        if record is None:
            outside = [item for item in selected if not self.assembly_record(item)]
            menu.addAction(f"Save selection as assembly… ({len(outside)} item{'s' if len(outside) != 1 else ''})",
                           lambda: self.main_window.save_selection_as_assembly(outside))
            return
        assembly_menu = menu.addMenu(f"Assembly “{record['name']}” ({'locked' if record['locked'] else 'unlocked'})")
        if record['locked']:
            assembly_menu.addAction("Unlock (edit, move or rotate parts)",
                                    lambda: self.set_assembly_locked(assembly_identifier, False))
        else:
            assembly_menu.addAction("Lock (move and edit as one block)",
                                    lambda: self.set_assembly_locked(assembly_identifier, True))
            outsiders = [item for item in selected if item.assembly_identifier != assembly_identifier]
            if outsiders:
                assembly_menu.addAction(f"Add {len(outsiders)} selected item{'s' if len(outsiders) != 1 else ''} "
                                        "to this assembly",
                                        lambda: self.add_items_to_assembly(assembly_identifier, outsiders))
            assembly_menu.addAction("Remove this item from the assembly",
                                    lambda: self.remove_items_from_assembly([clicked_item]))
        assembly_menu.addAction("Rename…", lambda: self.main_window.rename_assembly(assembly_identifier))
        assembly_menu.addAction("Save to library…", lambda: self.main_window.save_assembly(assembly_identifier))
        assembly_menu.addSeparator()
        assembly_menu.addAction("Dissolve (keep the parts)", lambda: self.dissolve_assembly(assembly_identifier))

    def block_menu(self, block, event):
        if not block.isSelected():
            self.clearSelection()
            block.setSelected(True)
        menu = QMenu()
        menu.addAction("Properties…", lambda: self.main_window.edit_block(block))
        if has_toggle(block.block_type):
            state_menu = menu.addMenu("State")
            for state_index, state_label in enumerate(block.type_information['states']):
                action = state_menu.addAction(state_label)
                action.setCheckable(True)
                action.setChecked(state_index == block.state)
                action.triggered.connect(lambda checked=False, state_index=state_index: block.set_state(state_index))
        rotate_targets = [item for item in self.selectedItems() if isinstance(item, BlockItem)] or [block]
        menu.addAction("Rotate right 90°\tR", lambda: [self.rotate_block(target, True) for target in rotate_targets])
        menu.addAction("Rotate left 90°\tShift+R",
                       lambda: [self.rotate_block(target, False) for target in rotate_targets])
        self._add_data_label_action(menu, block)
        if block.block_type == 'pump':
            menu.addAction("Pump curves…", lambda: self.main_window.edit_pump_curves(block))
            menu.addAction("Find lowest RPM…", lambda: self.main_window.find_lowest_speed(block.identifier))
        if block.block_type == 'orifice':
            menu.addAction("Size this orifice…", lambda: self.main_window.size_orifice(block.identifier))
        menu.addSeparator()
        selected_blocks = [item for item in self.selectedItems() if isinstance(item, BlockItem)]
        selected_items = self.selectedItems() or [block]
        menu.addAction("Copy\tCtrl+C", lambda: self.main_window.copy_items(selected_items))
        menu.addAction("Duplicate Block", lambda: self.duplicate_items(selected_items))
        menu.addAction("Save Block…", lambda: self.main_window.save_block(block))
        self._add_assembly_actions(menu, block)
        menu.addSeparator()
        if self.is_locked_item(block):
            menu.addAction("Delete assembly", lambda: self.delete_items(self.selectedItems()))
        else:
            menu.addAction("Delete Block", lambda: self.delete_items(selected_blocks or [block]))
        menu.exec(event.screenPos())

    def line_menu(self, line, event):
        control_held = bool(event.modifiers() & Qt.KeyboardModifier.ControlModifier)
        if not line.isSelected():
            if not control_held:
                self.clearSelection()
            line.setSelected(True)
        selected_lines = [item for item in self.selectedItems() if isinstance(item, LineItem)]
        menu = QMenu()
        if self.is_locked_item(line):
            # tubes inside a locked assembly can't be reshaped or removed on their own
            menu.addAction("Copy\tCtrl+C", lambda: self.main_window.copy_items(self.selectedItems() or [line]))
            menu.addAction("Properties…", lambda: self.main_window.edit_line(line))
            self._add_data_label_action(menu, line)
            self._add_assembly_actions(menu, line)
            menu.addSeparator()
            menu.addAction("Delete assembly", lambda: self.delete_items(self.selectedItems()))
            menu.exec(event.screenPos())
            return
        if len(selected_lines) >= 2:
            menu.addAction(f"Combine {len(selected_lines)} lines into one…",
                           lambda: self.combine_lines(selected_lines))
            menu.addAction("Auto-route selected lines", lambda: self.auto_route(selected_lines))
            menu.addAction("Delete selected lines", lambda: self.delete_items(selected_lines))
            menu.addSeparator()
        menu.addAction("Copy\tCtrl+C", lambda: self.main_window.copy_items(self.selectedItems() or [line]))
        menu.addAction("Properties…", lambda: self.main_window.edit_line(line))
        self._add_data_label_action(menu, line)
        _, closest_point, segment_index, _ = project_on_polyline(line.points, to_tuple(event.scenePos()))
        menu.addAction("Add segment here", lambda: self.add_segment(line, segment_index, closest_point))
        menu.addAction("Delete segment", lambda: self.delete_line_segment(line, segment_index))
        menu.addAction("Auto-route", lambda: self.auto_route([line]))
        menu.addAction("Reverse direction", lambda: self.reverse_line(line))
        self._add_assembly_actions(menu, line)
        menu.addSeparator()
        menu.addAction("Delete line", lambda: self.delete_items([line]))
        menu.exec(event.screenPos())

    def canvas_menu(self, event):
        menu = QMenu()
        add_menu = menu.addMenu("Add block here")
        position = to_tuple(event.scenePos())
        for block_type in BLOCK_ORDER:
            add_menu.addAction(BLOCK_TYPES[block_type]['label'],
                               lambda checked=False, block_type=block_type: self.add_block(
                                   {'type': block_type, 'properties': default_block_properties(block_type)},
                                   position))
        menu.addAction(self.main_window.labels_action)
        paste_action = menu.addAction("Paste here\tCtrl+V", lambda: self.main_window.paste_items(position))
        paste_action.setEnabled(self.main_window.clipboard_payload() is not None)
        menu.addAction("Auto-route all lines", lambda: self.auto_route(list(self.lines.values())))
        menu.addAction("Solve now", self.main_window.solve_now)
        menu.exec(event.screenPos())

    # ---- mouse in draw mode ----------------------------------------------------
    def mousePressEvent(self, event):
        if self.mode == 'draw':
            if event.button() == Qt.MouseButton.LeftButton:
                self.drawing_click(to_tuple(event.scenePos()))
            elif event.button() == Qt.MouseButton.RightButton:
                self.finish_drawing_with_free_end()
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self.mode == 'draw':
            self.drawing_move(to_tuple(event.scenePos()))
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.mode == 'draw':
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event):
        if self.mode == 'draw':
            self.finish_drawing_with_free_end()
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    # ---- drag & drop from the palette -----------------------------------------
    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat(BLOCK_MIME_TYPE):
            event.acceptProposedAction()
        else:
            super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        if event.mimeData().hasFormat(BLOCK_MIME_TYPE):
            event.acceptProposedAction()
        else:
            super().dragMoveEvent(event)

    def place_assembly(self, template, position):
        """Add a saved assembly from the library, centred at `position`, as a new locked assembly."""
        return self.paste_payload(template, position)

    def dropEvent(self, event):
        if event.mimeData().hasFormat(BLOCK_MIME_TYPE):
            template = json.loads(bytes(event.mimeData().data(BLOCK_MIME_TYPE).data()).decode('utf-8'))
            if template.get('format') == 'pylon_assembly':
                self.place_assembly(template, to_tuple(event.scenePos()))
            else:
                self.add_block(template_to_block_data(template), to_tuple(event.scenePos()))
            event.acceptProposedAction()
        else:
            super().dropEvent(event)

    # ---- model input / output -------------------------------------------------------
    def to_model(self):
        return {'format': 'pylon_schematic', 'version': FILE_FORMAT_VERSION,
                'application': f"{APPLICATION_NAME} {APPLICATION_VERSION}",
                'fluid': dict(self.fluid), 'settings': dict(self.settings),
                'blocks': [block.to_dictionary() for block in self.blocks.values()],
                'lines': [line.to_dictionary() for line in self.lines.values()],
                'scenarios': copy.deepcopy(self.scenarios),
                'active_scenario': self.active_scenario,
                'assemblies': [dict(record, identifier=identifier) for identifier, record in self.assemblies.items()],
                'view': {'show_data_labels': self.show_labels}}

    def clear_all(self):
        for line in list(self.lines.values()):
            line.remove_labels()
            self.removeItem(line)
        for block in list(self.blocks.values()):
            self.removeItem(block)
        for outline in self.assembly_outlines.values():
            self.removeItem(outline)
        self.blocks, self.lines = {}, {}
        self.assemblies, self.assembly_outlines = {}, {}
        self.results = None

    def load_model(self, model):
        model = upgrade_model(copy.deepcopy(model))
        for block_data in model.get('blocks', []):
            block_data.setdefault('x', 0)
            block_data.setdefault('y', 0)
        route_model_lines(model)
        self.loading = True
        try:
            self.clear_all()
            self.fluid = dict(DEFAULT_FLUID)
            self.fluid.update(model.get('fluid') or {})
            self.settings = dict(DEFAULT_SETTINGS)
            self.settings.update(model.get('settings') or {})
            self.scenarios = copy.deepcopy(model.get('scenarios') or [])
            self.active_scenario = model.get('active_scenario')
            self.show_labels = bool((model.get('view') or {}).get('show_data_labels', True))
            self.assemblies = {record['identifier']: {'name': record.get('name') or 'Assembly',
                                                      'locked': bool(record.get('locked', True))}
                               for record in (model.get('assemblies') or []) if record.get('identifier')}
            for block_data in model.get('blocks', []):
                if block_data.get('type') in BLOCK_TYPES:
                    self.add_block(block_data, select=False)
            for line_data in model.get('lines', []):
                if line_data.get('points') and len(line_data['points']) >= 2:
                    self.add_line(line_data)
        finally:
            self.loading = False
        for line in self.lines.values():
            port_moved = False
            for which_end in ('start', 'end'):
                attachment = getattr(line, which_end)
                if attachment and attachment.get('kind') == 'port' and attachment.get('block') not in self.blocks:
                    setattr(line, which_end, None)
                if attachment and attachment.get('kind') == 'line' and attachment.get('line') not in self.lines:
                    setattr(line, which_end, None)
                attachment = getattr(line, which_end)
                if attachment and attachment.get('kind') == 'port':
                    port_point, _ = self.attachment_point(attachment)
                    current = line.points[0] if which_end == 'start' else line.points[-1]
                    if port_point is not None and not points_close(port_point, current, 1.0):
                        port_moved = True
            if port_moved:
                # the port is no longer where the file says (e.g. a tank outlet saved before
                # it moved to the bottom of the tank): re-route the line to the port
                self.reroute_line(line)
            line.update_path()
        self.refresh_port_markers()
        self.notify()

    def apply_results(self, results):
        self.results = results
        for identifier, block in self.blocks.items():
            block.set_result(results['blocks'].get(identifier) if results else None)
        for identifier, line in self.lines.items():
            line.set_result(results['lines'].get(identifier) if results else None)

    def set_show_labels(self, visible):
        """Show or hide every flow / pressure label; clears the per-element overrides."""
        self.show_labels = bool(visible)
        for item in list(self.blocks.values()) + list(self.lines.values()):
            item.data_label_override = None
            item.apply_data_label_visibility()

    def apply_snapshot(self, snapshot):
        for identifier, saved in snapshot.get('blocks', {}).items():
            block = self.blocks.get(identifier)
            if block:
                block.state = normalized_state(block.block_type, saved.get('state', 0))
                block.properties.update(copy.deepcopy(saved.get('properties', {})))
                block.refresh()
        for identifier, saved in snapshot.get('lines', {}).items():
            line = self.lines.get(identifier)
            if line:
                line.properties.update(copy.deepcopy(saved.get('properties', {})))
                line.update_path()
        self.active_scenario = snapshot.get('name')
        self.notify()


def template_to_block_data(template):
    """Palette template / saved-block JSON -> block data for add_block."""
    template = upgrade_block_template(template) if 'props' in template else template
    block_type = template['type']
    properties = default_block_properties(block_type)
    properties.update(copy.deepcopy(template.get('properties') or {}))
    return {'type': block_type, 'properties': properties, 'state': int(template.get('state', 0)),
            'pump_curves': copy.deepcopy(template.get('pump_curves') or {}),
            'orientation': int(template.get('orientation', 0))}


# -----------------------------------------------------------------------------
# View
# -----------------------------------------------------------------------------
class SchematicView(QGraphicsView):
    def __init__(self, scene):
        super().__init__(scene)
        self.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing)
        self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        self.setAcceptDrops(True)
        self.setMouseTracking(True)
        self._pan_position = None

    def wheelEvent(self, event):
        wheel_steps = event.angleDelta().y() / 120.0
        zoom_factor = 1.15 ** wheel_steps
        current_zoom = self.transform().m11()
        if 0.15 < current_zoom * zoom_factor < 6:
            self.scale(zoom_factor, zoom_factor)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.MiddleButton:
            self._pan_position = event.position()
            self.viewport().setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._pan_position is not None:
            movement = event.position() - self._pan_position
            self._pan_position = event.position()
            self.horizontalScrollBar().setValue(int(self.horizontalScrollBar().value() - movement.x()))
            self.verticalScrollBar().setValue(int(self.verticalScrollBar().value() - movement.y()))
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.MiddleButton and self._pan_position is not None:
            self._pan_position = None
            scene = self.scene()
            mode = scene.mode if isinstance(scene, SchematicScene) else 'select'
            self.viewport().setCursor(Qt.CursorShape.ArrowCursor if mode == 'select'
                                      else Qt.CursorShape.CrossCursor)
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event):
        scene = self.scene()
        focus_item = scene.focusItem() if scene is not None else None
        if event.key() == Qt.Key.Key_Escape and isinstance(scene, SchematicScene) and focus_item is None:
            scene.escape()
            return
        super().keyPressEvent(event)

    def fit_all(self):
        scene = self.scene()
        if scene is None:
            return
        bounds = QRectF()
        for item in scene.items():
            if isinstance(item, (BlockItem, LineItem)):
                bounds = bounds.united(item.sceneBoundingRect())
        if bounds.isNull():
            return
        self.fitInView(bounds.adjusted(-60, -60, 60, 60), Qt.AspectRatioMode.KeepAspectRatio)
        if self.transform().m11() > 1.5:
            self.resetTransform()
            self.scale(1.5, 1.5)
            self.centerOn(bounds.center())


# -----------------------------------------------------------------------------
# Component palette (left column)
# -----------------------------------------------------------------------------
class PaletteEntry(QFrame):
    clicked = Signal(dict)
    delete_requested = Signal(str)

    def __init__(self, template, saved_path=None):
        super().__init__()
        self.template = template
        self.saved_path = saved_path
        self._press_position = None
        self.setObjectName('entry')
        self.setStyleSheet("#entry { background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
                           "stop:0 #14243c, stop:1 #070d18); border: 2px solid %(gold_middle)s; "
                           "border-radius: 10px; }"
                           "#entry:hover { border: 2px solid %(cyan)s; background: #10304c; }" % THEME)
        entry_layout = QHBoxLayout(self)
        entry_layout.setContentsMargins(4, 3, 6, 3)
        is_assembly = template.get('format') == 'pylon_assembly'
        icon = QLabel()
        icon.setPixmap(assembly_pixmap(template) if is_assembly
                       else symbol_pixmap(template['type'], int(template.get('state', 0))))
        icon.setFixedSize(78, 52)
        entry_layout.addWidget(icon)
        type_label = (f"Assembly · {len(template.get('blocks', []))} parts" if is_assembly
                      else BLOCK_TYPES[template['type']]['label'])
        text = QLabel(f"<b style='color:{THEME['heading']}'>{template.get('name') or type_label}</b>"
                      + (f"<br>{muted_html(type_label)}" if saved_path else ''))
        text.setWordWrap(True)
        entry_layout.addWidget(text, 1)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip("Click to add to the schematic, or drag it onto the canvas"
                        + ("\nRight-click to remove from the library" if saved_path else ''))

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._press_position = event.position()

    def mouseMoveEvent(self, event):
        if self._press_position is None or not (event.buttons() & Qt.MouseButton.LeftButton):
            return
        if (event.position() - self._press_position).manhattanLength() < QApplication.startDragDistance():
            return
        self._press_position = None
        drag = QDrag(self)
        mime_data = QMimeData()
        mime_data.setData(BLOCK_MIME_TYPE, json.dumps(self.template).encode('utf-8'))
        drag.setMimeData(mime_data)
        drag.setPixmap(assembly_pixmap(self.template) if self.template.get('format') == 'pylon_assembly'
                       else symbol_pixmap(self.template['type']))
        drag.exec(Qt.DropAction.CopyAction)

    def mouseReleaseEvent(self, event):
        if self._press_position is not None and event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit(self.template)
        self._press_position = None

    def contextMenuEvent(self, event):
        if not self.saved_path:
            return
        menu = QMenu(self)
        menu.addAction("Delete saved assembly" if self.template.get('format') == 'pylon_assembly'
                       else "Delete saved block", lambda: self.delete_requested.emit(self.saved_path))
        menu.exec(event.globalPos())


class PaletteWidget(QWidget):
    add_requested = Signal(dict)

    def __init__(self):
        super().__init__()
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        frame = BezelFrame("Component Blocks", padding=4)
        outer_layout.addWidget(frame)
        frame_layout = frame.content_layout
        self.scroll_area = QScrollArea()
        self.scroll_area.setObjectName('paletteScroll')
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll_area.setStyleSheet("#paletteScroll { border: none; background: transparent; }")
        self.inner_widget = QWidget()
        self.inner_widget.setObjectName('paletteInner')
        self.inner_widget.setStyleSheet("#paletteInner { background: transparent; }")
        self.entries_layout = QVBoxLayout(self.inner_widget)
        self.entries_layout.setContentsMargins(6, 6, 6, 6)
        self.entries_layout.setSpacing(5)
        self.scroll_area.setWidget(self.inner_widget)
        frame_layout.addWidget(self.scroll_area, 1)
        hint = QLabel("Click to add · drag onto canvas")
        hint.setStyleSheet(f"color:{THEME['muted']};")
        frame_layout.addWidget(hint)
        self.setMinimumWidth(240)
        self.refresh()

    def refresh(self):
        while self.entries_layout.count():
            layout_item = self.entries_layout.takeAt(0)
            widget = layout_item.widget()
            if widget is not None:
                widget.deleteLater()
        self.entries_layout.addWidget(self._header("Starting blocks"))
        for block_type in BLOCK_ORDER:
            entry = PaletteEntry({'type': block_type, 'properties': default_block_properties(block_type),
                                  'state': 0, 'name': BLOCK_TYPES[block_type]['label']})
            entry.clicked.connect(self.add_requested.emit)
            self.entries_layout.addWidget(entry)
        self.entries_layout.addWidget(self._header("Saved blocks & assemblies"))
        saved_blocks = list_saved_blocks()
        if not saved_blocks:
            placeholder = QLabel("<i>None yet — right-click a block → Save Block, or select connected "
                                 "parts → right-click → Save selection as assembly</i>")
            placeholder.setWordWrap(True)
            placeholder.setStyleSheet(f"color:{THEME['muted']};")
            self.entries_layout.addWidget(placeholder)
        for path, template in saved_blocks:
            entry = PaletteEntry(template, saved_path=path)
            entry.clicked.connect(self.add_requested.emit)
            entry.delete_requested.connect(self._delete_saved_block)
            self.entries_layout.addWidget(entry)
        self.entries_layout.addStretch(1)

    @staticmethod
    def _header(text):
        header = QLabel(text.upper())
        header.setStyleSheet(f"color:{THEME['heading']}; font-weight:600; font-size:8pt; padding-top:6px; "
                             f"letter-spacing:1px; border-bottom: 1px solid {THEME['gold_dark']};")
        return header

    def _delete_saved_block(self, path):
        if QMessageBox.question(self, "Delete saved block",
                                f"Delete {os.path.basename(path)} from the Component Blocks folder?"
                                ) == QMessageBox.StandardButton.Yes:
            try:
                os.remove(path)
            except OSError as error:
                QMessageBox.warning(self, "Delete saved block", str(error))
            self.refresh()


def list_saved_blocks():
    """(path, template) for every saved block, upgraded to the current keys."""
    saved_blocks = []
    if not os.path.isdir(BLOCK_DIRECTORY):
        return saved_blocks
    for file_name in sorted(os.listdir(BLOCK_DIRECTORY), key=str.lower):
        if not file_name.lower().endswith('.json'):
            continue
        path = os.path.join(BLOCK_DIRECTORY, file_name)
        try:
            with open(path, encoding='utf-8') as file_handle:
                template = json.load(file_handle)
            if template.get('format') == 'pylon_assembly' and isinstance(template.get('blocks'), list):
                saved_blocks.append((path, template))
            elif template.get('type') in BLOCK_TYPES:
                saved_blocks.append((path, upgrade_block_template(template)))
        except (OSError, ValueError):
            continue
    return saved_blocks


# -----------------------------------------------------------------------------
# Dialogs
# -----------------------------------------------------------------------------
class PropertyDialog(BezelDialog):
    def __init__(self, parent, title, definitions, values, calculated_values=None, name=None, states=None,
                 state=0, block_type=None, extra_buttons=None, tube_helper=False, note=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(500)
        self.definitions = {definition.key: definition for definition in definitions}
        self.editors = {}
        self.values = {}
        self.name = name
        self.state = state
        outer_layout = QVBoxLayout(self)
        if note:
            note_label = QLabel(note)
            note_label.setWordWrap(True)
            note_label.setStyleSheet(f"color:{THEME['muted']};")
            outer_layout.addWidget(note_label)
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        inner_widget = QWidget()
        form = QFormLayout(inner_widget)
        scroll_area.setWidget(inner_widget)
        outer_layout.addWidget(scroll_area, 1)
        self.name_editor = None
        if name is not None:
            self.name_editor = QLineEdit(name)
            form.addRow("Name", self.name_editor)
        self.state_selector = None
        if states:
            self.state_selector = QComboBox()
            self.state_selector.addItems(states)
            self.state_selector.setCurrentIndex(state)
            form.addRow("State (toggle)", self.state_selector)
        if tube_helper:
            tube_selector = QComboBox()
            tube_selector.addItem("— choose to fill outer / inner diameter —")
            tube_selector.addItems(TUBE_LABELS)
            tube_selector.currentIndexChanged.connect(
                lambda index: self._fill_tube(tube_selector.currentText()) if index > 0 else None)
            form.addRow("Standard tube size", tube_selector)
        calculated_header_added = False
        for definition in definitions:
            if definition.calculated:
                if calculated_values is None:
                    continue
                if not calculated_header_added:
                    form.addRow(QLabel(f"<b style='color:{THEME['heading']}'>CALCULATED</b>"))
                    calculated_header_added = True
                value = calculated_values.get(definition.key)
                display = QLineEdit(format_number(float(value), 6) if isinstance(value, (int, float)) else '—')
                display.setReadOnly(True)
                display.setStyleSheet(f"background:#0d1626; color:{THEME['cyan_soft']}; border-color:#2a3a4c;")
                form.addRow(self._row_label(definition), display)
                continue
            value = values.get(definition.key, definition.default)
            if definition.kind == 'tube':
                editor = QComboBox()
                editor.setEditable(True)
                editor.addItems(TUBE_LABELS)
                editor.setCurrentText(str(value))
            elif definition.kind == 'choice':
                editor = QComboBox()
                editor.addItems(definition.choices or [])
                editor.setCurrentText(str(value))
            elif definition.kind == 'int':
                editor = QSpinBox()
                editor.setRange(0, 9999)
                editor.setValue(int(round(float(value or 0))))
            else:
                editor = QLineEdit(f"{float(value):.6g}" if value is not None else '')
            if definition.tooltip:
                editor.setToolTip(definition.tooltip)
            form.addRow(self._row_label(definition), editor)
            self.editors[definition.key] = editor
        if block_type:
            for diameter_key, flow_coefficient_key, discharge_key in \
                    BLOCK_TYPES[block_type]['flow_coefficient_links']:
                if all(key in self.editors for key in (diameter_key, flow_coefficient_key, discharge_key)):
                    self.editors[diameter_key].textEdited.connect(
                        lambda _text, diameter=diameter_key, coefficient=flow_coefficient_key,
                        discharge=discharge_key: self._link(diameter, coefficient, discharge, False))
                    self.editors[discharge_key].textEdited.connect(
                        lambda _text, diameter=diameter_key, coefficient=flow_coefficient_key,
                        discharge=discharge_key: self._link(diameter, coefficient, discharge, False))
                    self.editors[flow_coefficient_key].textEdited.connect(
                        lambda _text, diameter=diameter_key, coefficient=flow_coefficient_key,
                        discharge=discharge_key: self._link(diameter, coefficient, discharge, True))
        if extra_buttons:
            button_row = QHBoxLayout()
            for text, handler in extra_buttons:
                button = QPushButton(text)
                button.clicked.connect(handler)
                button_row.addWidget(button)
            button_row.addStretch(1)
            outer_layout.addLayout(button_row)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        accent_buttons(button_box)
        button_box.accepted.connect(self._accept)
        button_box.rejected.connect(self.reject)
        outer_layout.addWidget(button_box)
        self.resize(540, min(760, 140 + 34 * (len(definitions) + 3)))

    @staticmethod
    def _row_label(definition):
        return f"{definition.label} [{definition.unit}]" if definition.unit else definition.label

    def _fill_tube(self, label):
        outer_diameter, wall_thickness = tube_outer_diameter_and_wall(label)
        if 'outer_diameter' in self.editors and 'inner_diameter' in self.editors:
            self.editors['outer_diameter'].setText(f"{outer_diameter:.6g}")
            self.editors['inner_diameter'].setText(f"{outer_diameter - 2 * wall_thickness:.6g}")

    def _link(self, diameter_key, flow_coefficient_key, discharge_key, from_flow_coefficient):
        """Keep diameter and Cv consistent while typing (Cv = 29.84 Cd d^2, [ISA])."""
        try:
            discharge_coefficient = float(self.editors[discharge_key].text())
            if from_flow_coefficient:
                self.editors[diameter_key].setText(f"{orifice_diameter_from_flow_coefficient(float(self.editors[flow_coefficient_key].text()), discharge_coefficient):.6g}")
            else:
                self.editors[flow_coefficient_key].setText(f"{flow_coefficient_from_orifice_diameter(float(self.editors[diameter_key].text()), discharge_coefficient):.6g}")
        except ValueError:
            pass

    def _accept(self):
        collected = {}
        for key, editor in self.editors.items():
            definition = self.definitions[key]
            if isinstance(editor, QComboBox):
                collected[key] = editor.currentText().strip()
                if definition.kind == 'tube' and parse_tube(collected[key]) is None:
                    QMessageBox.warning(self, "Invalid tube size",
                                        f"{definition.label}: use a standard size or 'outer diameter x wall' "
                                        f'in inches, e.g. 5/16" x 0.035"')
                    editor.setFocus()
                    return
            elif isinstance(editor, QSpinBox):
                collected[key] = editor.value()
            elif isinstance(editor, QLineEdit):
                try:
                    collected[key] = float(editor.text())
                except ValueError:
                    QMessageBox.warning(self, "Invalid number",
                                        f"{definition.label}: '{editor.text()}' is not a number.")
                    editor.setFocus()
                    return
        if ('inner_diameter' in collected and 'outer_diameter' in collected
                and not (0 < collected['inner_diameter'] < collected['outer_diameter'])):
            QMessageBox.warning(self, "Invalid tube", "Inner diameter must be > 0 and smaller than the outer diameter.")
            return
        self.values = collected
        if self.name_editor is not None:
            self.name = self.name_editor.text().strip() or self.name
        if self.state_selector is not None:
            self.state = self.state_selector.currentIndex()
        self.accept()


# ordinal blue ramp for the dark surface: low speed darker, high speed brighter
ORDINAL_BLUES = ['#256abf', '#2a78d6', '#5598e7', '#86b6ef', '#b7d3f6', '#cde2fb']


class CurvePlot(QWidget):
    """Pump curves: pressure rise on x, flow on y, one line per speed (ordinal blue ramp)."""

    def __init__(self):
        super().__init__()
        self.pump_curves = {}
        self.selected_speed = None
        self.setMinimumSize(420, 280)
        self.setMouseTracking(True)
        self._plotted_points = []

    def set_curves(self, pump_curves, selected_speed=None):
        self.pump_curves = pump_curves
        self.selected_speed = selected_speed
        self.update()

    @staticmethod
    def _colors(count):
        if count <= 1:
            return [ORDINAL_BLUES[2]]
        return [ORDINAL_BLUES[round(index * (len(ORDINAL_BLUES) - 1) / (count - 1))] for index in range(count)]

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(THEME['background_field']))
        self._plotted_points = []
        if not self.pump_curves:
            painter.setPen(COLOR_MUTED)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No curves yet — use Import CSV…")
            return
        speeds = sorted(self.pump_curves, key=float)
        pressure_values = [point[0] for speed in speeds for point in self.pump_curves[speed]] + [0.0]
        flow_values = [point[1] for speed in speeds for point in self.pump_curves[speed]] + [0.0]
        maximum_pressure = max(pressure_values) * 1.05 or 1.0
        maximum_flow = max(flow_values) * 1.08 or 1.0
        left_margin, right_margin, top_margin, bottom_margin = 62, 90, 14, 44
        plot_width = self.width() - left_margin - right_margin
        plot_height = self.height() - top_margin - bottom_margin

        def to_x(value):
            return left_margin + value / maximum_pressure * plot_width

        def to_y(value):
            return top_margin + plot_height - value / maximum_flow * plot_height
        painter.setFont(make_font(7.5))
        for tick_index in range(6):
            pressure_tick, flow_tick = maximum_pressure * tick_index / 5, maximum_flow * tick_index / 5
            painter.setPen(QPen(QColor('#1b2a3d'), 1))
            painter.drawLine(QPointF(to_x(pressure_tick), top_margin), QPointF(to_x(pressure_tick), top_margin + plot_height))
            painter.drawLine(QPointF(left_margin, to_y(flow_tick)), QPointF(left_margin + plot_width, to_y(flow_tick)))
            painter.setPen(COLOR_MUTED)
            painter.drawText(QRectF(to_x(pressure_tick) - 30, top_margin + plot_height + 4, 60, 14),
                             Qt.AlignmentFlag.AlignCenter, format_number(pressure_tick, 3))
            painter.drawText(QRectF(4, to_y(flow_tick) - 7, left_margin - 8, 14), Qt.AlignmentFlag.AlignRight,
                             format_number(flow_tick, 3))
        painter.setPen(QPen(QColor(THEME['gold']), 1.2))
        painter.drawLine(QPointF(left_margin, top_margin + plot_height),
                         QPointF(left_margin + plot_width, top_margin + plot_height))
        painter.setPen(COLOR_TEXT)
        painter.drawText(QRectF(left_margin, top_margin + plot_height + 22, plot_width, 16),
                         Qt.AlignmentFlag.AlignCenter, "Pump ΔP [psid]")
        painter.save()
        painter.translate(14, top_margin + plot_height / 2)
        painter.rotate(-90)
        painter.drawText(QRectF(-plot_height / 2, -8, plot_height, 16), Qt.AlignmentFlag.AlignCenter, "Flow [mL/min]")
        painter.restore()
        for speed, color in zip(speeds, self._colors(len(speeds))):
            points = sorted(self.pump_curves[speed])
            painter.setPen(QPen(QColor(color), 3.0 if speed == self.selected_speed else 2.0))
            painter.drawPolyline(QPolygonF([QPointF(to_x(pressure), to_y(flow)) for pressure, flow in points]))
            painter.setBrush(QColor(color))
            for pressure, flow in points:
                painter.drawEllipse(QPointF(to_x(pressure), to_y(flow)), 3, 3)
                self._plotted_points.append((to_x(pressure), to_y(flow), speed, pressure, flow))
            last_pressure, last_flow = points[-1]
            painter.setPen(COLOR_TEXT)
            painter.drawText(QRectF(to_x(last_pressure) + 6, to_y(last_flow) - 7, right_margin - 8, 14),
                             Qt.AlignmentFlag.AlignLeft, f"{speed} rpm")

    def mouseMoveEvent(self, event):
        cursor = event.position()
        nearest = None
        for x, y, speed, pressure, flow in self._plotted_points:
            distance = math.hypot(x - cursor.x(), y - cursor.y())
            if distance < 12 and (nearest is None or distance < nearest[0]):
                nearest = (distance, speed, pressure, flow)
        if nearest:
            QToolTip.showText(event.globalPosition().toPoint(),
                              f"{nearest[1]} rpm\nΔP {format_number(nearest[2])} psid\n"
                              f"Flow {format_number(nearest[3])} mL/min", self)
        else:
            QToolTip.hideText()


class PumpCurveDialog(BezelDialog):
    def __init__(self, parent, block, density):
        super().__init__(parent)
        self.block = block
        self.density = density
        self.pump_curves = copy.deepcopy(block.pump_curves)
        self.setWindowTitle(f"Pump curves — {block.name}")
        self.resize(900, 520)
        dialog_layout = QHBoxLayout(self)
        left_column = QVBoxLayout()
        left_column.addWidget(QLabel("<b>RPM curves</b>"))
        self.speed_list = QListWidget()
        self.speed_list.currentTextChanged.connect(self._select)
        left_column.addWidget(self.speed_list, 1)
        for text, handler in (("Import CSV…", self._import), ("Add / edit RPM…", self._add_speed),
                              ("Delete curve", self._delete_curve), ("Save CSV template…", self._save_template)):
            button = QPushButton(text)
            button.clicked.connect(handler)
            left_column.addWidget(button)
        dialog_layout.addLayout(left_column)
        middle_column = QVBoxLayout()
        middle_column.addWidget(QLabel("Points (ΔP psid, flow mL/min) — editable"))
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["ΔP [psid]", "Flow [mL/min]"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.itemChanged.connect(self._table_edited)
        middle_column.addWidget(self.table, 1)
        point_buttons = QHBoxLayout()
        add_button = QPushButton("Add point")
        add_button.clicked.connect(self._add_point)
        remove_button = QPushButton("Remove point")
        remove_button.clicked.connect(self._remove_point)
        point_buttons.addWidget(add_button)
        point_buttons.addWidget(remove_button)
        middle_column.addLayout(point_buttons)
        dialog_layout.addLayout(middle_column)
        right_column = QVBoxLayout()
        self.plot = CurvePlot()
        right_column.addWidget(self.plot, 1)
        note = QLabel("One curve per RPM. When the pump is ON it follows the curve for its set speed: flow vs. "
                      "ΔP is linear between points (extrapolated beyond the ends), and a speed between two "
                      "stored curves is interpolated between them.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{THEME['muted']};")
        right_column.addWidget(note)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        accent_buttons(button_box)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        right_column.addWidget(button_box)
        dialog_layout.addLayout(right_column, 1)
        self._loading = False
        self._refresh_list()

    def _refresh_list(self, selected_speed=None):
        self.speed_list.blockSignals(True)
        self.speed_list.clear()
        speeds = sorted(self.pump_curves, key=float)
        self.speed_list.addItems(speeds)
        self.speed_list.blockSignals(False)
        if speeds:
            speed = selected_speed if selected_speed in speeds else speeds[0]
            self.speed_list.setCurrentRow(speeds.index(speed))
            self._select(speed)
        else:
            self._select(None)

    def _select(self, speed):
        self._loading = True
        self.table.setRowCount(0)
        if speed and speed in self.pump_curves:
            for pressure, flow in self.pump_curves[speed]:
                row = self.table.rowCount()
                self.table.insertRow(row)
                self.table.setItem(row, 0, QTableWidgetItem(f"{pressure:.6g}"))
                self.table.setItem(row, 1, QTableWidgetItem(f"{flow:.6g}"))
        self._loading = False
        self.plot.set_curves(self.pump_curves, speed)

    def _current_speed(self):
        item = self.speed_list.currentItem()
        return item.text() if item else None

    def _table_edited(self, _item):
        if self._loading:
            return
        speed = self._current_speed()
        if not speed:
            return
        points = []
        for row in range(self.table.rowCount()):
            try:
                points.append([float(self.table.item(row, 0).text()), float(self.table.item(row, 1).text())])
            except (AttributeError, ValueError):
                continue
        self.pump_curves[speed] = sorted(points)
        self.plot.set_curves(self.pump_curves, speed)

    def _add_point(self):
        speed = self._current_speed()
        if not speed:
            return
        points = self.pump_curves.get(speed, [])
        last_point = points[-1] if points else [0.0, 0.0]
        points.append([last_point[0] + 50.0, last_point[1]])
        self.pump_curves[speed] = sorted(points)
        self._select(speed)

    def _remove_point(self):
        speed = self._current_speed()
        row = self.table.currentRow()
        if speed and 0 <= row < len(self.pump_curves.get(speed, [])) and len(self.pump_curves[speed]) > 2:
            del self.pump_curves[speed][row]
            self._select(speed)

    def _add_speed(self):
        value, accepted = QInputDialog.getDouble(self, "Add RPM curve", "RPM:", 3000, 1, 1e6, 0)
        if accepted:
            speed = str(int(round(value)))
            self.pump_curves.setdefault(speed, [[0.0, 0.0], [100.0, 0.0]])
            self._refresh_list(speed)

    def _delete_curve(self):
        speed = self._current_speed()
        if speed:
            self.pump_curves.pop(speed, None)
            self._refresh_list()

    def _import(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Import pump curve CSV", SCRIPT_DIRECTORY,
                                                "CSV files (*.csv *.txt);;All files (*)")
        last_imported_speed = None
        for path in paths:
            def ask_speed():
                value, accepted = QInputDialog.getDouble(self, "RPM for curve",
                                                         f"{os.path.basename(path)} has no 'rpm' column.\nRPM:",
                                                         3000, 1, 1e6, 0)
                return value if accepted else None
            try:
                imported = parse_pump_curve_csv(path, self.density, ask_speed)
            except (OSError, ValueError) as error:
                QMessageBox.warning(self, "Import failed", f"{os.path.basename(path)}:\n{error}")
                continue
            self.pump_curves.update(imported)
            last_imported_speed = sorted(imported, key=float)[-1]
        if last_imported_speed:
            self._refresh_list(last_imported_speed)

    def _save_template(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save pump-curve template",
                                              os.path.join(SCRIPT_DIRECTORY, 'pump_curve_template.csv'),
                                              "CSV files (*.csv)")
        if not path:
            return
        with open(path, 'w', newline='', encoding='utf-8') as file_handle:
            writer = csv.writer(file_handle)
            writer.writerow(['rpm', 'dp_psid', 'flow_ml_min'])
            for speed in (2000, 3000):
                for pressure_rise in (0, 100, 200, 300, 400):
                    writer.writerow([speed, pressure_rise, round(0.12 * speed - 0.12 * pressure_rise, 2)])

    def accept(self):
        too_short = [speed for speed, points in self.pump_curves.items() if len(points) < 2]
        if too_short:
            QMessageBox.warning(self, "Pump curves",
                                f"Curves need at least two points: {', '.join(too_short)} rpm")
            return
        super().accept()


class SettingsDialog(BezelDialog):
    def __init__(self, parent, fluid, settings):
        super().__init__(parent)
        self.setWindowTitle("Fluid and loss coefficients")
        self.resize(560, 640)
        dialog_layout = QVBoxLayout(self)
        fluid_form = QFormLayout()
        dialog_layout.addWidget(QLabel("<b>Fluid</b> (incompressible, adiabatic)"))
        self.preset_selector = QComboBox()
        self.preset_selector.addItems(list(FLUID_PRESETS) + ['Custom'])
        self.preset_selector.setCurrentText(fluid.get('name') if fluid.get('name') in FLUID_PRESETS else 'Custom')
        self.density_editor = QLineEdit(f"{float(fluid['density']):.6g}")
        self.viscosity_editor = QLineEdit(f"{float(fluid['viscosity_millipascal_seconds']):.6g}")
        self.preset_selector.currentTextChanged.connect(self._apply_preset)
        fluid_form.addRow("Preset", self.preset_selector)
        fluid_form.addRow("Density [kg/m³]", self.density_editor)
        fluid_form.addRow("Dynamic viscosity [mPa·s]", self.viscosity_editor)
        dialog_layout.addLayout(fluid_form)
        note = QLabel("Presets are approximate values near 20 °C — enter your own for other temperatures.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{THEME['muted']};")
        dialog_layout.addWidget(note)
        dialog_layout.addWidget(QLabel("<b>Loss coefficients and analysis settings</b>"))
        settings_form = QFormLayout()
        self.setting_editors = {}
        for key, label, default, tooltip in SETTING_DEFINITIONS:
            editor = QLineEdit(f"{float(settings.get(key, default)):.6g}")
            editor.setToolTip(f"{tooltip}\nDefault: {default}")
            settings_form.addRow(label, editor)
            self.setting_editors[key] = editor
        dialog_layout.addLayout(settings_form)
        restore_button = QPushButton("Restore defaults")
        restore_button.clicked.connect(self._restore_defaults)
        dialog_layout.addWidget(restore_button)
        dialog_layout.addStretch(1)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        accent_buttons(button_box)
        button_box.accepted.connect(self._accept)
        button_box.rejected.connect(self.reject)
        dialog_layout.addWidget(button_box)
        self.fluid, self.settings = None, None

    def _restore_defaults(self):
        for key, _, default, _ in SETTING_DEFINITIONS:
            self.setting_editors[key].setText(f"{default:.6g}")

    def _apply_preset(self, name):
        if name in FLUID_PRESETS:
            density, viscosity = FLUID_PRESETS[name]
            self.density_editor.setText(f"{density:.6g}")
            self.viscosity_editor.setText(f"{viscosity:.6g}")

    def _accept(self):
        try:
            density = float(self.density_editor.text())
            viscosity = float(self.viscosity_editor.text())
            settings = {key: float(editor.text()) for key, editor in self.setting_editors.items()}
        except ValueError:
            QMessageBox.warning(self, "Invalid number", "All fields must be numbers.")
            return
        if density <= 0 or viscosity <= 0:
            QMessageBox.warning(self, "Invalid fluid", "Density and viscosity must be positive.")
            return
        name = self.preset_selector.currentText()
        if name in FLUID_PRESETS and FLUID_PRESETS[name] != (density, viscosity):
            name = f"{name} (edited)"
        self.fluid = {'name': name, 'density': density, 'viscosity_millipascal_seconds': viscosity}
        self.settings = settings
        self.accept()


def make_table(headers):
    table = QTableWidget(0, len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    table.horizontalHeader().setStretchLastSection(True)
    table.verticalHeader().setVisible(False)
    return table


def add_table_row(table, cells, background_color=None):
    row = table.rowCount()
    table.insertRow(row)
    for column, value in enumerate(cells):
        item = QTableWidgetItem(value if isinstance(value, str) else format_number(value))
        if background_color:
            item.setBackground(QColor(background_color))
        table.setItem(row, column, item)


class ScenarioPicker(QWidget):
    """'Current state' plus a checklist of saved scenarios."""

    def __init__(self, scenarios):
        super().__init__()
        picker_layout = QVBoxLayout(self)
        picker_layout.setContentsMargins(0, 0, 0, 0)
        self.scenario_list = QListWidget()
        current_item = QListWidgetItem("Current state (as shown on the canvas)")
        current_item.setFlags(current_item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        current_item.setCheckState(Qt.CheckState.Checked if not scenarios else Qt.CheckState.Unchecked)
        self.scenario_list.addItem(current_item)
        self.snapshots = [None]
        for scenario in scenarios:
            scenario_item = QListWidgetItem(f"Scenario: {scenario['name']}")
            scenario_item.setFlags(scenario_item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            scenario_item.setCheckState(Qt.CheckState.Checked)
            self.scenario_list.addItem(scenario_item)
            self.snapshots.append(scenario)
        self.scenario_list.setMaximumHeight(130)
        picker_layout.addWidget(self.scenario_list)

    def chosen_snapshots(self, current_model):
        chosen = []
        for index in range(self.scenario_list.count()):
            if self.scenario_list.item(index).checkState() == Qt.CheckState.Checked:
                snapshot = self.snapshots[index]
                chosen.append(snapshot if snapshot is not None else make_snapshot(current_model, 'Current state'))
        return chosen


def block_name(scene, identifier):
    return scene.blocks[identifier].name if identifier in scene.blocks else identifier


class OrificeSizingDialog(BezelDialog):
    def __init__(self, main_window, orifice_identifier=None):
        super().__init__(main_window)
        self.main_window = main_window
        self.setWindowTitle("Bypass orifice sizing")
        self.resize(1000, 600)
        scene = main_window.scene
        dialog_layout = QVBoxLayout(self)
        introduction = QLabel(
            "Finds the diameter of an orifice placed in parallel with the pump (one side on the "
            "pump-outlet manifold, the other on the tank-outlet / pump-inlet manifold) so that "
            "<b>pump flow = thruster flow + orifice back flow</b> with every ON thruster at its "
            "nominal flow (its flow at nominal inlet pressure).")
        introduction.setWordWrap(True)
        dialog_layout.addWidget(introduction)
        form = QFormLayout()
        self.orifice_selector = QComboBox()
        self.orifice_identifiers = [block.identifier for block in sorted(scene.blocks.values(), key=lambda candidate: candidate.name)
                                    if block.block_type == 'orifice']
        self.orifice_selector.addItems([scene.blocks[identifier].name for identifier in self.orifice_identifiers])
        if orifice_identifier in self.orifice_identifiers:
            self.orifice_selector.setCurrentIndex(self.orifice_identifiers.index(orifice_identifier))
        form.addRow("Orifice block", self.orifice_selector)
        self.scenario_picker = ScenarioPicker(scene.scenarios)
        form.addRow("Evaluate", self.scenario_picker)
        self.tolerance_editor = QLineEdit(f"{scene.settings.get('acceptance_tolerance_percent', 2.0):.4g}")
        form.addRow("Thruster-flow tolerance [%]", self.tolerance_editor)
        dialog_layout.addLayout(form)
        button_row = QHBoxLayout()
        calculate_button = QPushButton("Calculate")
        calculate_button.clicked.connect(self.run)
        self.apply_button = QPushButton("Apply recommended diameter")
        self.apply_button.setEnabled(False)
        self.apply_button.clicked.connect(self.apply)
        button_row.addWidget(calculate_button)
        button_row.addWidget(self.apply_button)
        button_row.addStretch(1)
        dialog_layout.addLayout(button_row)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        dialog_layout.addWidget(self.summary)
        self.table = make_table(["Scenario", "Ideal diameter [in]", "Thruster flow at recommended diameter [g/s]",
                                 "Nominal demand [g/s]", "Error [%]", "Pump flow [g/s]",
                                 "Orifice back flow [g/s]", "Pump ΔP [psid]", "Pump outlet [psia]",
                                 "Limit warnings"])
        dialog_layout.addWidget(self.table, 1)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        accent_buttons(button_box)
        button_box.rejected.connect(self.reject)
        dialog_layout.addWidget(button_box)
        self.best_diameter = None
        if not self.orifice_identifiers:
            self.summary.setText(warning_html("Add an Orifice block to the schematic first."))
            calculate_button.setEnabled(False)

    def run(self):
        orifice_identifier = self.orifice_identifiers[self.orifice_selector.currentIndex()]
        model = self.main_window.scene.to_model()
        snapshots = self.scenario_picker.chosen_snapshots(model)
        if not snapshots:
            return
        try:
            tolerance_percent = float(self.tolerance_editor.text())
        except ValueError:
            tolerance_percent = 2.0
        progress_dialog = QProgressDialog("Sizing orifice…", "Cancel", 0, 100, self)
        progress_dialog.setWindowModality(Qt.WindowModality.WindowModal)
        progress_dialog.setMinimumDuration(200)

        class Cancelled(Exception):
            pass

        def report_progress(percent):
            progress_dialog.setValue(percent)
            QApplication.processEvents()
            if progress_dialog.wasCanceled():
                raise Cancelled()
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            sizing = size_orifice_for_all_scenarios(model, orifice_identifier, snapshots, tolerance_percent,
                                                    progress=report_progress)
        except Cancelled:
            self.summary.setText("Cancelled.")
            return
        finally:
            QApplication.restoreOverrideCursor()
            progress_dialog.close()
        self.table.setRowCount(0)
        if sizing['best_diameter'] is None:
            self.summary.setText(warning_html(sizing['message']))
            return
        self.best_diameter = sizing['best_diameter']
        self.apply_button.setEnabled(True)
        worst_error_percent = max(abs(row['evaluation']['relative_error']) for row in sizing['rows']) * 100
        window = sizing['window']
        window_text = (f"Diameters from <b>{window[0]:.5f}</b> to <b>{window[1]:.5f} in</b> keep every scenario "
                       f"within ±{tolerance_percent:g}% of nominal thruster flow with no limit warnings."
                       if window else warning_html(f"No single diameter keeps every scenario within "
                                                   f"±{tolerance_percent:g}% without limit warnings — see the table."))
        self.summary.setText(f"<b>Recommended orifice diameter: {self.best_diameter:.5f} in</b> "
                             f"({self.best_diameter * 25.4:.4f} mm) — worst-case thruster-flow error "
                             f"{worst_error_percent:.2f}%.<br>{window_text}")
        scene = self.main_window.scene
        for row in sizing['rows']:
            evaluation = row['evaluation']
            warnings_text = '; '.join(f"{block_name(scene, identifier)}: {text}"
                                      for identifier, text in evaluation['violations'])
            acceptable = abs(evaluation['relative_error']) * 100 <= tolerance_percent and not evaluation['violations']
            add_table_row(self.table, [row['name'],
                                       f"{row['ideal_diameter']:.5f}" if row['ideal_diameter'] else row['ideal_message'],
                                       evaluation['thruster_flow'] / KILOGRAMS_PER_GRAM,
                                       evaluation['demand'] / KILOGRAMS_PER_GRAM,
                                       evaluation['relative_error'] * 100,
                                       evaluation['pump_flow'] / KILOGRAMS_PER_GRAM,
                                       evaluation['orifice_flow'] / KILOGRAMS_PER_GRAM,
                                       evaluation['pump_pressure_rise'], evaluation['pump_outlet_pressure'],
                                       warnings_text or '—'], None if acceptable else THEME['row_warning'])

    def apply(self):
        orifice_identifier = self.orifice_identifiers[self.orifice_selector.currentIndex()]
        scene = self.main_window.scene
        block = scene.blocks.get(orifice_identifier)
        if block and self.best_diameter:
            block.properties['orifice_diameter'] = round(self.best_diameter, 6)
            for scenario in scene.scenarios:
                if orifice_identifier in scenario.get('blocks', {}):
                    scenario['blocks'][orifice_identifier]['properties']['orifice_diameter'] = round(self.best_diameter, 6)
            block.refresh()
            scene.notify()
            self.main_window.scenario_panel.refresh()
            self.summary.setText(self.summary.text() + "<br><i>Applied to the orifice and all saved scenarios.</i>")


class LowestSpeedDialog(BezelDialog):
    def __init__(self, main_window, pump_identifier=None):
        super().__init__(main_window)
        self.main_window = main_window
        self.setWindowTitle("Lowest pump RPM")
        self.resize(980, 560)
        scene = main_window.scene
        dialog_layout = QVBoxLayout(self)
        introduction = QLabel("Runs each scenario with the pump ON at every imported "
                              "RPM (ascending) and reports the lowest RPM at which every ON thruster receives at "
                              "least its nominal flow (minus the tolerance).")
        introduction.setWordWrap(True)
        dialog_layout.addWidget(introduction)
        form = QFormLayout()
        self.pump_selector = QComboBox()
        self.pump_identifiers = [block.identifier for block in sorted(scene.blocks.values(), key=lambda candidate: candidate.name)
                                 if block.block_type == 'pump']
        self.pump_selector.addItems([f"{scene.blocks[identifier].name} "
                                     f"({len(scene.blocks[identifier].pump_curves)} curves)"
                                     for identifier in self.pump_identifiers])
        if pump_identifier in self.pump_identifiers:
            self.pump_selector.setCurrentIndex(self.pump_identifiers.index(pump_identifier))
        form.addRow("Pump", self.pump_selector)
        self.scenario_picker = ScenarioPicker(scene.scenarios)
        form.addRow("Evaluate", self.scenario_picker)
        self.tolerance_editor = QLineEdit(f"{scene.settings.get('acceptance_tolerance_percent', 2.0):.4g}")
        form.addRow("Thruster-flow tolerance [%]", self.tolerance_editor)
        dialog_layout.addLayout(form)
        calculate_button = QPushButton("Calculate")
        calculate_button.clicked.connect(self.run)
        dialog_layout.addWidget(calculate_button)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        dialog_layout.addWidget(self.summary)
        self.table = make_table(["Scenario", "RPM", "Pump flow [g/s]", "Thruster flow [g/s]",
                                 "Nominal demand [g/s]", "Minimum thruster / nominal", "Pump ΔP [psid]",
                                 "Meets demand", "Limit warnings"])
        dialog_layout.addWidget(self.table, 1)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        accent_buttons(button_box)
        button_box.rejected.connect(self.reject)
        dialog_layout.addWidget(button_box)
        if not self.pump_identifiers:
            self.summary.setText(warning_html("Add a Pump block first."))
            calculate_button.setEnabled(False)

    def run(self):
        pump_identifier = self.pump_identifiers[self.pump_selector.currentIndex()]
        scene = self.main_window.scene
        if not scene.blocks[pump_identifier].pump_curves:
            self.summary.setText(warning_html("This pump has no curves — right-click it → Pump curves… → "
                                              "Import CSV."))
            return
        model = scene.to_model()
        snapshots = self.scenario_picker.chosen_snapshots(model)
        try:
            tolerance_percent = float(self.tolerance_editor.text())
        except ValueError:
            tolerance_percent = 2.0
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            scenario_results = find_lowest_pump_speed(model, pump_identifier, snapshots, tolerance_percent)
        finally:
            QApplication.restoreOverrideCursor()
        summary_lines = []
        self.table.setRowCount(0)
        for scenario_result in scenario_results:
            name = scenario_result['name']
            if scenario_result['demand'] <= 0:
                summary_lines.append(f"<b>{name}</b>: no thruster is ON.")
            elif scenario_result['lowest'] is None:
                summary_lines.append(f"<b>{name}</b>: " + warning_html(
                    f"no curve RPM meets the demand ({scenario_result['demand'] / KILOGRAMS_PER_GRAM:.3g} g/s)."))
            else:
                extra = ''
                if scenario_result['lowest_without_warnings'] is None:
                    extra = " — " + warning_html("but pressure limits are exceeded at every RPM that meets it")
                elif scenario_result['lowest_without_warnings'] != scenario_result['lowest']:
                    extra = (f" (lowest with no limit warnings: "
                             f"{scenario_result['lowest_without_warnings']:.0f} rpm)")
                summary_lines.append(f"<b>{name}</b>: lowest RPM <b>{scenario_result['lowest']:.0f}</b>{extra}")
            for row in scenario_result['rows']:
                warnings_text = '; '.join(f"{block_name(scene, identifier)}: {text}"
                                          for identifier, text in row['violations'])
                row_color = (THEME['row_good'] if row['speed'] == scenario_result['lowest']
                             else (None if row['meets_demand'] else THEME['row_dimmed']))
                add_table_row(self.table, [name, f"{row['speed']:.0f}", row['pump_flow'] / KILOGRAMS_PER_GRAM,
                                           row['thruster_flow'] / KILOGRAMS_PER_GRAM,
                                           row['demand'] / KILOGRAMS_PER_GRAM, row['minimum_flow_ratio'],
                                           row['pump_pressure_rise'], 'yes' if row['meets_demand'] else 'no',
                                           warnings_text or '—'], row_color)
        self.summary.setText('<br>'.join(summary_lines))


# -----------------------------------------------------------------------------
# Side panels
# -----------------------------------------------------------------------------
class ResultsPanel(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        panel_layout = QVBoxLayout(self)
        panel_layout.setContentsMargins(0, 4, 0, 0)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        panel_layout.addWidget(BezelFrame("System status", self.summary))
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.table = make_table(["Item", "Type", "State", "P in [psia]", "P out [psia]", "ΔP [psid]",
                                 "ṁ [g/s]", "Notes"])
        self.table.setAlternatingRowColors(True)
        self.table.cellDoubleClicked.connect(self._go_to_row)
        splitter.addWidget(BezelFrame("Component data", self.table))
        self.warning_list = QListWidget()
        self.warning_list.itemDoubleClicked.connect(
            lambda item: self._select_item(item.data(Qt.ItemDataRole.UserRole)))
        splitter.addWidget(BezelFrame("Warnings", self.warning_list))
        splitter.setSizes([420, 150])
        panel_layout.addWidget(splitter, 1)
        self._row_identifiers = []

    def _select_item(self, identifier):
        scene = self.main_window.scene
        item = scene.blocks.get(identifier) or scene.lines.get(identifier)
        if item is not None:
            scene.clearSelection()
            item.setSelected(True)
            for view in scene.views():
                view.centerOn(item)

    def _go_to_row(self, row, _column):
        if 0 <= row < len(self._row_identifiers):
            self._select_item(self._row_identifiers[row])

    def show_results(self, results):
        scene = self.main_window.scene
        self.table.setRowCount(0)
        self._row_identifiers = []
        if results is None:
            self.summary.setText("")
            return
        density = results['density']
        thrusters = [block for block in scene.blocks.values() if block.block_type == 'thruster' and block.state == 0]
        pumps = [block for block in scene.blocks.values() if block.block_type == 'pump' and block.state == 0]
        orifices = [block for block in scene.blocks.values() if block.block_type == 'orifice' and block.state == 0]
        thruster_flow = sum(results['blocks'][block.identifier].get('mass_flow', 0.0) for block in thrusters)
        nominal_flow = sum(thruster_nominal_mass_flow(block.properties, density) for block in thrusters)
        total_thrust = sum(results['blocks'][block.identifier]['calculated'].get('calculated_thrust', 0.0)
                           for block in thrusters)
        status_color = THEME['success'] if results['converged'] else THEME['warning']
        summary_lines = [f"<span style='color:{status_color}'><b>{results['message']}</b></span>"
                         f" · fluid {scene.fluid['name']} (ρ {format_number(scene.fluid['density'])} kg/m³)"]
        if scene.active_scenario:
            summary_lines.append(f"Scenario: <b>{scene.active_scenario}</b>")
        if thrusters:
            deviation_percent = (thruster_flow / nominal_flow - 1) * 100 if nominal_flow else 0
            summary_lines.append(f"Thrusters ON: {len(thrusters)} · flow "
                                 f"<b>{format_number(thruster_flow / KILOGRAMS_PER_GRAM, 4)} g/s</b> "
                                 f"(nominal {format_number(nominal_flow / KILOGRAMS_PER_GRAM, 4)} g/s, "
                                 f"{format_number(deviation_percent, 3)}%) · thrust "
                                 f"<b>{format_number(total_thrust, 4)} N</b>")
        if pumps:
            pump_flow = sum(results['blocks'][block.identifier].get('mass_flow', 0.0) for block in pumps)
            orifice_flow = sum(results['blocks'][block.identifier].get('mass_flow', 0.0) for block in orifices)
            summary_lines.append(f"Pump flow {format_number(pump_flow / KILOGRAMS_PER_GRAM, 4)} g/s = thruster "
                                 f"{format_number(thruster_flow / KILOGRAMS_PER_GRAM, 4)} + orifice "
                                 f"{format_number(orifice_flow / KILOGRAMS_PER_GRAM, 4)} g/s (balance "
                                 f"{format_number((pump_flow - thruster_flow - orifice_flow) / KILOGRAMS_PER_GRAM, 3)})")
        self.summary.setText('<br>'.join(summary_lines))
        for block in sorted(scene.blocks.values(), key=lambda item: item.name):
            block_result = results['blocks'].get(block.identifier) or {}
            ports = block_result.get('ports', {})
            calculated = block_result.get('calculated', {})
            inlet_psia = ports.get('inlet', ports.get('outlet', float('nan'))) / PASCALS_PER_PSI
            outlet_psia = (ports.get('outlet', float('nan')) / PASCALS_PER_PSI
                           if 'outlet' in ports and 'inlet' in ports else float('nan'))
            pressure_drop = inlet_psia - outlet_psia if not math.isnan(outlet_psia) else float('nan')
            mass_flow = block_result.get('mass_flow')
            note = block_result.get('label', '')
            if block.block_type == 'transducer':
                note = f"reads {format_number(calculated.get('calculated_pressure'), 5)} psia"
            if block.block_type == 'pump':
                pressure_drop = -calculated.get('calculated_pressure_rise', float('nan'))
            if block_result.get('warnings'):
                note = '⚠ ' + '; '.join(block_result['warnings'])
            add_table_row(self.table, [block.name, block.type_information['label'],
                                       block.type_information['states'][block.state]
                                       if has_toggle(block.block_type) else '—', inlet_psia, outlet_psia,
                                       pressure_drop,
                                       (mass_flow / KILOGRAMS_PER_GRAM) if mass_flow is not None else '—', note],
                          THEME['row_warning'] if block_result.get('warnings') else None)
            self._row_identifiers.append(block.identifier)
        for line in sorted(scene.lines.values(), key=lambda item: item.name):
            line_result = results['lines'].get(line.identifier) or {}
            for section_index, section in enumerate(line_result.get('sections', [])):
                row_name = line.name if len(line_result['sections']) == 1 else f"{line.name}.{section_index + 1}"
                add_table_row(self.table, [row_name, 'Line', '', section['start_pressure'] / PASCALS_PER_PSI,
                                           section['end_pressure'] / PASCALS_PER_PSI,
                                           section['pressure_drop'] / PASCALS_PER_PSI,
                                           section['mass_flow'] / KILOGRAMS_PER_GRAM,
                                           f"velocity {format_number(section['velocity'], 3)} m/s · Reynolds "
                                           f"{format_number(section['reynolds_number'], 3)}"])
                self._row_identifiers.append(line.identifier)
        self.warning_list.clear()
        for identifier, text in results['warnings']:
            item = QListWidgetItem(f"⚠ {block_name(scene, identifier)}: {text}")
            item.setData(Qt.ItemDataRole.UserRole, identifier)
            item.setForeground(COLOR_WARNING)
            self.warning_list.addItem(item)
        if not results['warnings']:
            self.warning_list.addItem("No warnings.")


class ScenarioPanel(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 4, 0, 0)
        frame = BezelFrame("Scenarios")
        outer_layout.addWidget(frame)
        frame_layout = frame.content_layout
        information = QLabel("A scenario stores every block's toggle state and properties, and every "
                             "line's properties. Double-click to apply.")
        information.setWordWrap(True)
        information.setStyleSheet(f"color:{THEME['muted']};")
        frame_layout.addWidget(information)
        self.scenario_list = QListWidget()
        self.scenario_list.itemDoubleClicked.connect(lambda _: self.apply_selected())
        frame_layout.addWidget(self.scenario_list, 1)
        for text, handler in (("Save current as new scenario…", self.save_new),
                              ("Update selected from current", self.update_selected),
                              ("Apply selected", self.apply_selected), ("Rename…", self.rename_selected),
                              ("Duplicate", self.duplicate_selected), ("Delete", self.delete_selected)):
            button = QPushButton(text)
            button.clicked.connect(handler)
            frame_layout.addWidget(button)

    def refresh(self):
        scene = self.main_window.scene
        self.scenario_list.clear()
        for scenario in scene.scenarios:
            item = QListWidgetItem(scenario['name'])
            if scenario['name'] == scene.active_scenario:
                font = item.font()
                font.setBold(True)
                item.setFont(font)
                item.setText(f"{scenario['name']}   ◀ active")
            item.setData(Qt.ItemDataRole.UserRole, scenario['name'])
            self.scenario_list.addItem(item)

    def _selected_scenario(self):
        item = self.scenario_list.currentItem()
        if not item:
            return None
        name = item.data(Qt.ItemDataRole.UserRole)
        for scenario in self.main_window.scene.scenarios:
            if scenario['name'] == name:
                return scenario
        return None

    def _unique_name(self, base_name):
        names = {scenario['name'] for scenario in self.main_window.scene.scenarios}
        if base_name not in names:
            return base_name
        counter = 2
        while f"{base_name} ({counter})" in names:
            counter += 1
        return f"{base_name} ({counter})"

    def save_new(self):
        scene = self.main_window.scene
        name, accepted = QInputDialog.getText(self, "New scenario", "Scenario name:",
                                              text=self._unique_name(f"Scenario {len(scene.scenarios) + 1}"))
        if accepted and name.strip():
            snapshot = make_snapshot(scene.to_model(), self._unique_name(name.strip()))
            scene.scenarios.append(snapshot)
            scene.active_scenario = snapshot['name']
            self.refresh()
            self.main_window.set_modified()

    def update_selected(self):
        scenario = self._selected_scenario()
        if scenario:
            scene = self.main_window.scene
            replacement = make_snapshot(scene.to_model(), scenario['name'])
            scenario.clear()
            scenario.update(replacement)
            scene.active_scenario = scenario['name']
            self.refresh()
            self.main_window.set_modified()

    def apply_selected(self):
        scenario = self._selected_scenario()
        if scenario:
            self.main_window.scene.apply_snapshot(scenario)
            self.refresh()

    def rename_selected(self):
        scenario = self._selected_scenario()
        if scenario:
            name, accepted = QInputDialog.getText(self, "Rename scenario", "Name:", text=scenario['name'])
            if accepted and name.strip() and name.strip() != scenario['name']:
                was_active = self.main_window.scene.active_scenario == scenario['name']
                scenario['name'] = self._unique_name(name.strip())
                if was_active:
                    self.main_window.scene.active_scenario = scenario['name']
                self.refresh()
                self.main_window.set_modified()

    def duplicate_selected(self):
        scenario = self._selected_scenario()
        if scenario:
            duplicate = copy.deepcopy(scenario)
            duplicate['name'] = self._unique_name(scenario['name'] + ' copy')
            self.main_window.scene.scenarios.append(duplicate)
            self.refresh()
            self.main_window.set_modified()

    def delete_selected(self):
        scenario = self._selected_scenario()
        if scenario and QMessageBox.question(self, "Delete scenario", f"Delete scenario '{scenario['name']}'?") \
                == QMessageBox.StandardButton.Yes:
            scene = self.main_window.scene
            scene.scenarios.remove(scenario)
            if scene.active_scenario == scenario['name']:
                scene.active_scenario = None
            self.refresh()
            self.main_window.set_modified()


# -----------------------------------------------------------------------------
# Main window
# -----------------------------------------------------------------------------
HELP_TEXT = """
<h3>Pylon</h3>
<p><b>Ports</b> — blocks take flow in on the left and out on the right (tee branch and 90° fitting
outlet at the bottom; the cross has branches at the top and bottom). The tank has its inlet on top and
its outlet on the bottom.</p>
<p><b>Blocks</b> — click or drag from the <i>Component Blocks</i> column. Drag a block to move it.
Click the <b>centre</b> of a block to toggle its state (tees, crosses and 90° fittings have no toggle).
Right-click for Properties, State, Rotate, Copy, Duplicate Block, Save Block and Delete Block.
Double-click the body for properties.
<b>Assemblies</b> — select connected parts and tubes, right-click → <i>Save selection as assembly…</i>
to save them to the Component Blocks library as one block. A <b>locked</b> assembly (padlock closed)
selects, moves, copies and deletes as one unit and protects its parts; right-click → Assembly →
<i>Unlock</i> to move, rotate or edit parts, add or remove parts, then <i>Lock</i> and <i>Save to library</i>.
<b>Ctrl+Z</b> undoes and <b>Ctrl+Y</b> redoes. <b>R</b> / <b>Shift+R</b> (or the toolbar's Rotate
right / Rotate left) turn the selected blocks 90°.
<b>D</b> (toolbar: Show flow / pressure data) shows or hides every flow and pressure label; while
they are hidden, hover over a block or line to see its data, and right-click an item →
<i>Show flow / pressure data here</i> to show or hide that item's labels on their own.
<b>Ctrl+C</b> copies the selected blocks together with the lines between them; <b>Ctrl+V</b>
(or Ctrl+P) pastes a copy next to the original, or right-click the canvas → Paste here.
Ctrl+D duplicates the selection in one step.
An ON pump follows its imported RPM curves; its speed (rpm) is shown on the block — click the number to
edit it. A pump with no curves is solved as stopped and flagged.</p>
<p><b>Lines</b> — drag from a port to another port (or onto a line) for an auto-routed tube that
avoids every component. Or press <b>L</b> (Draw Line) and click horizontal/vertical sections; click
a port or line to finish, right-click/double-click to finish with a free end, Esc to cancel.
Drag a segment to move it, drag an end to reconnect. Right-click a line to add or delete segments,
auto-route, reverse or edit properties. Select several lines (Ctrl-click or rubber band),
right-click → <b>Combine</b>.</p>
<p><b>Results</b> — mass flow and inlet→outlet ΔP are shown beside every tube (arrows show flow
direction); block results beside each block. Middle-drag pans, wheel zooms, F fits.</p>
<p><b>Analysis</b> — Scenarios tab stores toggle states and properties. Analysis → Size bypass
orifice… and Find lowest pump RPM… work on the current state and/or any saved scenarios.</p>
"""


class MainWindow(QMainWindow):
    def __init__(self, path=None):
        super().__init__()
        self.resize(1500, 900)
        self.file_path = None
        self.modified = False
        self._initial_guess = None
        self._clipboard = None
        self._paste_count = 0
        self._history = []
        self._history_index = -1
        self._restoring_history = False
        self._history_timer = QTimer(self)
        self._history_timer.setSingleShot(True)
        self._history_timer.setInterval(120)
        self._history_timer.timeout.connect(self._record_history)
        self.scene = SchematicScene(self)
        self.view = SchematicView(self.scene)
        self.setCentralWidget(BezelFrame("Schematic", self.view, rim=8.0, radius=16.0, padding=4))
        self.scene.model_changed.connect(self._on_model_changed)
        self.scene.message.connect(lambda text: self.statusBar().showMessage(text, 8000))
        self.solve_timer = QTimer(self)
        self.solve_timer.setSingleShot(True)
        self.solve_timer.setInterval(150)
        self.solve_timer.timeout.connect(self.solve_now)

        self.block_palette = PaletteWidget()
        self.block_palette.add_requested.connect(self.add_from_palette)
        palette_dock = QDockWidget("Components", self)
        palette_dock.setObjectName('paletteDock')
        palette_dock.setWidget(self.block_palette)
        palette_dock.setFeatures(QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        palette_dock.setTitleBarWidget(QWidget())
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, palette_dock)

        self.results_panel = ResultsPanel(self)
        self.scenario_panel = ScenarioPanel(self)
        results_dock = QDockWidget("Results", self)
        results_dock.setObjectName('resultsDock')
        results_dock.setWidget(self.results_panel)
        scenario_dock = QDockWidget("Scenarios", self)
        scenario_dock.setObjectName('scenarioDock')
        scenario_dock.setWidget(self.scenario_panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, results_dock)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, scenario_dock)
        for dock in (results_dock, scenario_dock):
            dock.setFeatures(QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
            dock.setTitleBarWidget(QWidget())
        self.setTabPosition(Qt.DockWidgetArea.RightDockWidgetArea, QTabWidget.TabPosition.North)
        self.tabifyDockWidget(results_dock, scenario_dock)
        results_dock.raise_()
        self.resizeDocks([results_dock], [520], Qt.Orientation.Horizontal)
        self._build_actions(results_dock, scenario_dock)
        if path:
            self.open_file(path)
        else:
            self.load_model(build_example_model(), None)
        QTimer.singleShot(0, self.view.fit_all)

    # ---- actions --------------------------------------------------------------
    def _make_action(self, text, handler, shortcut=None, checkable=False, tooltip=None):
        action = QAction(text, self)
        if shortcut:
            action.setShortcut(QKeySequence(shortcut))
        action.setCheckable(checkable)
        if tooltip:
            action.setStatusTip(tooltip)
            action.setToolTip(tooltip)
        action.triggered.connect(handler)
        return action

    def _build_actions(self, results_dock, scenario_dock):
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("&File")
        file_menu.addAction(self._make_action("New", self.new_file, "Ctrl+N"))
        file_menu.addAction(self._make_action("Open…", self.open_dialog, "Ctrl+O"))
        file_menu.addAction(self._make_action("Open example", lambda: self._confirm_discard() and
                                              self.load_model(build_example_model(), None)))
        file_menu.addAction(self._make_action("Save", self.save, "Ctrl+S"))
        file_menu.addAction(self._make_action("Save As…", self.save_as, "Ctrl+Shift+S"))
        file_menu.addSeparator()
        file_menu.addAction(self._make_action("Export image (PNG)…", self.export_image))
        file_menu.addAction(self._make_action("Open Component Blocks folder location", self.show_block_directory))
        file_menu.addSeparator()
        file_menu.addAction(self._make_action("Quit", self.close, "Ctrl+Q"))

        edit_menu = menu_bar.addMenu("&Edit")
        self.undo_action = self._make_action("Undo", self.undo, QKeySequence.StandardKey.Undo,
                                             tooltip="Undo the last change (Ctrl+Z)")
        self.redo_action = self._make_action("Redo", self.redo, None, tooltip="Redo (Ctrl+Y or Ctrl+Shift+Z)")
        self.redo_action.setShortcuts([QKeySequence("Ctrl+Y"), QKeySequence("Ctrl+Shift+Z")])
        edit_menu.addAction(self.undo_action)
        edit_menu.addAction(self.redo_action)
        edit_menu.addSeparator()
        self.select_action = self._make_action("Select / Move", lambda: self.set_mode('select'), "S", True,
                                               "Select, move and edit (drag from a port to connect)")
        self.draw_action = self._make_action("Draw Line", lambda: self.set_mode('draw'), "L", True,
                                             "Draw a plumbing line in horizontal/vertical sections")
        mode_group = QActionGroup(self)
        mode_group.addAction(self.select_action)
        mode_group.addAction(self.draw_action)
        self.select_action.setChecked(True)
        edit_menu.addAction(self.select_action)
        edit_menu.addAction(self.draw_action)
        edit_menu.addSeparator()
        edit_menu.addAction(self._make_action("Delete selected", self.delete_selected, QKeySequence.StandardKey.Delete))
        edit_menu.addAction(self._make_action("Copy", self.copy_selected, QKeySequence.StandardKey.Copy,
                                              tooltip="Copy the selected blocks and the lines between them"))
        paste_action = self._make_action("Paste", lambda: self.paste_items(None), None,
                                         tooltip="Paste copied blocks and lines (Ctrl+V or Ctrl+P)")
        paste_action.setShortcuts([QKeySequence(QKeySequence.StandardKey.Paste), QKeySequence("Ctrl+P")])
        edit_menu.addAction(paste_action)
        edit_menu.addAction(self._make_action("Duplicate selection", self.duplicate_selected, "Ctrl+D"))
        self.rotate_left_action = self._make_action("⟲ Rotate left", lambda: self.rotate_selected(False), "Shift+R",
                                                    tooltip="Rotate the selected blocks 90° counter-clockwise (Shift+R)")
        self.rotate_right_action = self._make_action("Rotate right ⟳", lambda: self.rotate_selected(True), "R",
                                                     tooltip="Rotate the selected blocks 90° clockwise (R)")
        edit_menu.addAction(self.rotate_left_action)
        edit_menu.addAction(self.rotate_right_action)
        edit_menu.addAction(self._make_action("Select all", self.select_all, "Ctrl+A"))

        view_menu = menu_bar.addMenu("&View")
        view_menu.addAction(self._make_action("Fit to window", self.view.fit_all, "F"))
        view_menu.addAction(self._make_action("Zoom in", lambda: self.view.scale(1.2, 1.2), "Ctrl+="))
        view_menu.addAction(self._make_action("Zoom out", lambda: self.view.scale(1 / 1.2, 1 / 1.2), "Ctrl+-"))
        self.labels_action = self._make_action("Show flow / pressure data", self.set_show_data_labels, "D", True,
                                               "Show or hide the flow and pressure labels on every block and "
                                               "line (D). When hidden, hover over an item to see its data; "
                                               "right-click an item to show or hide its own labels.")
        self.labels_action.setChecked(True)
        view_menu.addAction(self.labels_action)
        view_menu.addSeparator()
        view_menu.addAction(results_dock.toggleViewAction())
        view_menu.addAction(scenario_dock.toggleViewAction())

        analysis_menu = menu_bar.addMenu("&Analysis")
        analysis_menu.addAction(self._make_action("Solve now", self.solve_now, "F5"))
        self.auto_solve_action = self._make_action("Auto-solve on every change",
                                                   lambda enabled: enabled and self.solve_now(), None, True)
        self.auto_solve_action.setChecked(True)
        analysis_menu.addAction(self.auto_solve_action)
        analysis_menu.addSeparator()
        analysis_menu.addAction(self._make_action("Size bypass orifice…", lambda: self.size_orifice(None)))
        analysis_menu.addAction(self._make_action("Find lowest pump RPM…", lambda: self.find_lowest_speed(None)))
        analysis_menu.addSeparator()
        analysis_menu.addAction(self._make_action("Fluid and loss coefficients…", self.edit_settings))

        help_menu = menu_bar.addMenu("&Help")
        help_menu.addAction(self._make_action(
            "Quick guide", lambda: QMessageBox.information(self, "Quick guide", HELP_TEXT), "F1"))
        help_menu.addAction(self._make_action("About", lambda: QMessageBox.about(
            self, "About", f"<b>{APPLICATION_NAME}</b> {APPLICATION_VERSION}<br>Steady incompressible "
                           "flow-network solver for spacecraft propulsion feed systems.")))

        toolbar = self.addToolBar("Main")
        toolbar.setObjectName('mainToolbar')
        toolbar.setMovable(False)
        toolbar.addAction(self.undo_action)
        toolbar.addAction(self.redo_action)
        toolbar.addSeparator()
        toolbar.addAction(self.select_action)
        toolbar.addAction(self.draw_action)
        toolbar.addAction(self.rotate_left_action)
        toolbar.addAction(self.rotate_right_action)
        toolbar.addSeparator()
        toolbar.addAction(self.labels_action)
        toolbar.addSeparator()
        toolbar.addAction(self._make_action("Solve", self.solve_now))
        toolbar.addAction(self.auto_solve_action)
        toolbar.addAction(self._make_action("Fit", self.view.fit_all))
        toolbar.addSeparator()
        toolbar.addAction(self._make_action("Fluid…", self.edit_settings))
        toolbar.addAction(self._make_action("Size Orifice…", lambda: self.size_orifice(None)))
        toolbar.addAction(self._make_action("Lowest RPM…", lambda: self.find_lowest_speed(None)))

    def set_mode(self, mode):
        self.scene.set_mode(mode)
        (self.draw_action if mode == 'draw' else self.select_action).setChecked(True)

    def select_all(self):
        for item in self.scene.items():
            if isinstance(item, (BlockItem, LineItem)):
                item.setSelected(True)

    # ---- model / solve -----------------------------------------------------------
    def _on_model_changed(self, physics_changed):
        self.set_modified()
        if physics_changed and self.auto_solve_action.isChecked():
            self.solve_timer.start()
        self.scenario_panel.refresh()

    def solve_now(self):
        model = self.scene.to_model()
        try:
            results = solve_model(model, initial_guess=self._initial_guess)
            if not results['converged'] and self._initial_guess is not None:
                results = solve_model(model)
        except Exception as error:  # keep the editor alive whatever the network looks like
            import traceback
            traceback.print_exc()
            self.statusBar().showMessage(f"Solver error: {error}", 10000)
            return
        self._initial_guess = results['initial_guess'] if results['converged'] else None
        self.scene.apply_results(results)
        self.results_panel.show_results(results)
        warning_count = len(results['warnings'])
        self.statusBar().showMessage(f"{results['message']}"
                                     + (f" · {warning_count} warning(s)" if warning_count else ''), 6000)

    def load_model(self, model, path):
        self.scene.load_model(model)
        self.labels_action.setChecked(self.scene.show_labels)
        self.file_path = path
        self._initial_guess = None
        self.scenario_panel.refresh()
        self.solve_now()
        self.modified = False
        self._update_title()
        self.reset_history()
        QTimer.singleShot(0, self.view.fit_all)

    def set_modified(self, modified=True):
        self.modified = modified
        self._update_title()
        if modified:
            self._request_history_record()

    def _update_title(self):
        name = os.path.basename(self.file_path) if self.file_path else 'Untitled'
        self.setWindowTitle(f"{name}{' *' if self.modified else ''} — {APPLICATION_NAME}")

    # ---- palette ------------------------------------------------------------------
    def add_from_palette(self, template):
        view_center = self.view.mapToScene(self.view.viewport().rect().center())
        if template.get('format') == 'pylon_assembly':
            self.scene.place_assembly(template, (snap_to_grid(view_center.x()), snap_to_grid(view_center.y())))
            return
        width, height = block_size(template['type'])
        occupied = self.scene.block_rectangles()
        base_position = (snap_to_grid(view_center.x()), snap_to_grid(view_center.y()))
        chosen_position = base_position
        for attempt in range(400):
            offset_x = ((attempt % 7) - 3) * 160
            offset_y = ((attempt // 7) % 7 - 3) * 140
            candidate = (base_position[0] + offset_x, base_position[1] + offset_y) if attempt else base_position
            footprint = (candidate[0] - width / 2 - 30, candidate[1] - height / 2 - 40,
                         candidate[0] + width / 2 + 30, candidate[1] + height / 2 + 40)
            if not any(footprint[0] < other[2] and other[0] < footprint[2]
                       and footprint[1] < other[3] and other[1] < footprint[3] for other in occupied):
                chosen_position = candidate
                break
        self.scene.add_block(template_to_block_data(template), chosen_position)

    # ---- assemblies --------------------------------------------------------------
    def _write_library_file(self, name, template, title):
        safe_name = re.sub(r'[^\w\- .()]+', '_', name.strip()).strip() or 'assembly'
        os.makedirs(BLOCK_DIRECTORY, exist_ok=True)
        path = os.path.join(BLOCK_DIRECTORY, safe_name + '.json')
        if os.path.exists(path) and QMessageBox.question(
                self, title, f"'{safe_name}.json' already exists. Overwrite?") != QMessageBox.StandardButton.Yes:
            return None
        try:
            with open(path, 'w', encoding='utf-8') as file_handle:
                json.dump(template, file_handle, indent=2)
        except OSError as error:
            QMessageBox.warning(self, title, str(error))
            return None
        self.block_palette.refresh()
        self.statusBar().showMessage(f"Saved to {path}", 8000)
        return path

    def assembly_template(self, assembly_identifier, name):
        """Library template for an assembly: its parts and tubes, all in one locked assembly."""
        member_blocks, member_lines = self.scene.assembly_members(assembly_identifier)
        payload = self.scene.copy_payload(member_blocks + member_lines)
        for item in payload['blocks'] + payload['lines']:
            item['assembly'] = 'assembly'
        return {'format': 'pylon_assembly', 'version': FILE_FORMAT_VERSION, 'name': name,
                'blocks': payload['blocks'], 'lines': payload['lines'],
                'assemblies': {'assembly': {'name': name, 'locked': True}}}

    def save_selection_as_assembly(self, items):
        items = [item for item in items if isinstance(item, (BlockItem, LineItem))]
        if not any(isinstance(item, BlockItem) for item in items):
            self.statusBar().showMessage("Select the parts (and the lines between them) to save as an assembly.", 6000)
            return
        taken = {record['name'] for record in self.scene.assemblies.values()}
        default_name = incremented_name('Assembly', taken) if 'Assembly' in taken else 'Assembly'
        name, accepted = QInputDialog.getText(self, "Save assembly", "Name for the assembly:", text=default_name)
        if not accepted or not name.strip():
            return
        # connecting tubes between the selected parts belong to the assembly too
        payload = self.scene.copy_payload(items)
        line_identifiers = {line['identifier'] for line in payload['lines']} if payload else set()
        members = items + [line for line in self.scene.lines.values()
                           if line.identifier in line_identifiers and line not in items]
        assembly_identifier = self.scene.create_assembly(members, name.strip(), locked=True)
        self._write_library_file(name.strip(), self.assembly_template(assembly_identifier, name.strip()),
                                 "Save assembly")

    def save_assembly(self, assembly_identifier):
        record = self.scene.assemblies.get(assembly_identifier)
        if record is None:
            return
        name, accepted = QInputDialog.getText(self, "Save assembly to library", "Name for the saved assembly:",
                                              text=record['name'])
        if not accepted or not name.strip():
            return
        self.scene.rename_assembly(assembly_identifier, name.strip())
        self._write_library_file(name.strip(), self.assembly_template(assembly_identifier, name.strip()),
                                 "Save assembly")

    def rename_assembly(self, assembly_identifier):
        record = self.scene.assemblies.get(assembly_identifier)
        if record is None:
            return
        name, accepted = QInputDialog.getText(self, "Rename assembly", "Assembly name:", text=record['name'])
        if accepted and name.strip():
            self.scene.rename_assembly(assembly_identifier, name.strip())

    def _unlock_for_editing(self, item):
        """Parts of a locked assembly can't be modified; offer to unlock it."""
        record = self.scene.assembly_record(item)
        if not record or not record['locked']:
            return True
        if QMessageBox.question(self, "Locked assembly",
                                f"“{record['name']}” is locked. Unlock it to modify its parts?"
                                ) != QMessageBox.StandardButton.Yes:
            return False
        self.scene.set_assembly_locked(item.assembly_identifier, False)
        return True

    def save_block(self, block):
        name, accepted = QInputDialog.getText(self, "Save Block", "Name for the saved block:", text=block.name)
        if not accepted or not name.strip():
            return
        safe_name = re.sub(r'[^\w\- .()]+', '_', name.strip()).strip() or 'block'
        os.makedirs(BLOCK_DIRECTORY, exist_ok=True)
        path = os.path.join(BLOCK_DIRECTORY, safe_name + '.json')
        if os.path.exists(path) and QMessageBox.question(
                self, "Save Block", f"'{safe_name}.json' already exists. Overwrite?") != QMessageBox.StandardButton.Yes:
            return
        template = {'format': 'pylon_block', 'version': FILE_FORMAT_VERSION, 'name': name.strip(),
                    'type': block.block_type, 'state': block.state, 'orientation': block.orientation,
                    'properties': copy.deepcopy(block.properties), 'pump_curves': copy.deepcopy(block.pump_curves)}
        try:
            with open(path, 'w', encoding='utf-8') as file_handle:
                json.dump(template, file_handle, indent=2)
        except OSError as error:
            QMessageBox.warning(self, "Save Block", str(error))
            return
        self.block_palette.refresh()
        self.statusBar().showMessage(f"Saved block to {path}", 8000)

    def show_block_directory(self):
        os.makedirs(BLOCK_DIRECTORY, exist_ok=True)
        QMessageBox.information(self, "Component Blocks folder", BLOCK_DIRECTORY)

    # ---- editing -------------------------------------------------------------------
    def edit_block(self, block):
        if not self._unlock_for_editing(block):
            return
        calculated_values = (block.result or {}).get('calculated', {})
        extra_buttons = []
        if block.block_type == 'pump':
            extra_buttons.append(("Pump curves…", lambda: self.edit_pump_curves(block)))
        dialog = PropertyDialog(self, f"{block.type_information['label']} properties",
                                block.type_information['properties'], block.properties,
                                calculated_values=calculated_values, name=block.name,
                                states=block.type_information['states'] if has_toggle(block.block_type) else None,
                                state=block.state,
                                block_type=block.block_type, extra_buttons=extra_buttons)
        if dialog.exec():
            block.properties.update(dialog.values)
            update_linked_flow_coefficients(block.block_type, block.properties)
            block.name = dialog.name
            block.state = dialog.state
            block.refresh()
            self.scene.refresh_port_markers()
            self.scene.notify()

    def edit_line(self, line):
        if not self._unlock_for_editing(line):
            return
        calculated_values = (line.result or {}).get('calculated', {})
        dialog = PropertyDialog(self, "Plumbing line properties", LINE_PROPERTIES, line.properties,
                                calculated_values=calculated_values, name=line.name, tube_helper=True,
                                note="Bends and length are used in the flow calculation only — they are "
                                     "not tied to the drawn shape.")
        if dialog.exec():
            line.properties.update(dialog.values)
            line.name = dialog.name
            line.update_path()
            self.scene.notify()

    def edit_pump_curves(self, block):
        if not self._unlock_for_editing(block):
            return
        dialog = PumpCurveDialog(self, block, float(self.scene.fluid['density']))
        if dialog.exec():
            block.pump_curves = dialog.pump_curves
            stored_speeds = sorted(float(speed) for speed in block.pump_curves if block.pump_curves[speed])
            speed = float(block.properties.get('revolutions_per_minute', 0))
            if stored_speeds and not stored_speeds[0] <= speed <= stored_speeds[-1]:
                # keep the set speed inside the imported curves (interpolation range)
                block.properties['revolutions_per_minute'] = min(stored_speeds, key=lambda stored: abs(stored - speed))
                self.statusBar().showMessage(
                    f"{block.name}: speed set to {block.properties['revolutions_per_minute']:g} rpm "
                    "(inside the imported curves).", 8000)
            block.refresh()
            self.scene.notify()

    def edit_settings(self):
        dialog = SettingsDialog(self, self.scene.fluid, self.scene.settings)
        if dialog.exec():
            self.scene.fluid = dialog.fluid
            self.scene.settings.update(dialog.settings)
            self.scene.notify()

    def delete_selected(self):
        if self.scene.focusItem() is not None:
            return
        self.scene.delete_items(self.scene.selectedItems())

    def duplicate_selected(self):
        if self._text_editing():
            return
        selected_items = [item for item in self.scene.selectedItems() if isinstance(item, (BlockItem, LineItem))]
        if selected_items:
            self.scene.duplicate_items(selected_items)

    # ---- copy / paste --------------------------------------------------------------
    CLIPBOARD_MIME_TYPE = 'application/x-pylon-schematic'

    def _text_editing(self):
        """True while the keyboard belongs to a text field (so Ctrl+C/V/P stay with it)."""
        if self.scene.focusItem() is not None:
            return True
        focus_widget = QApplication.focusWidget()
        return isinstance(focus_widget, (QLineEdit, QAbstractSpinBox, QTextEdit, QPlainTextEdit))

    def copy_selected(self):
        if self._text_editing():
            return
        self.copy_items(self.scene.selectedItems())

    def copy_items(self, items):
        payload = self.scene.copy_payload([item for item in items if isinstance(item, (BlockItem, LineItem))])
        if payload is None:
            self.statusBar().showMessage("Nothing selected to copy.", 4000)
            return
        self._clipboard = payload
        self._paste_count = 0
        try:                                   # also on the system clipboard, for other windows
            mime_data = QMimeData()
            mime_data.setData(self.CLIPBOARD_MIME_TYPE, json.dumps(payload).encode('utf-8'))
            QApplication.clipboard().setMimeData(mime_data)
        except Exception:
            pass
        self.statusBar().showMessage(f"Copied {len(payload['blocks'])} block(s) and "
                                     f"{len(payload['lines'])} line(s).", 5000)

    def clipboard_payload(self):
        try:
            mime_data = QApplication.clipboard().mimeData()
            if mime_data is not None and mime_data.hasFormat(self.CLIPBOARD_MIME_TYPE):
                payload = json.loads(bytes(mime_data.data(self.CLIPBOARD_MIME_TYPE)).decode('utf-8'))
                if isinstance(payload, dict) and payload.get('format') in ('pylon_clipboard', 'propulsion_clipboard'):
                    if payload != self._clipboard:
                        self._clipboard, self._paste_count = payload, 0
                    return payload
        except Exception:
            pass
        return self._clipboard

    def paste_items(self, position=None):
        if position is None and self._text_editing():
            return
        payload = self.clipboard_payload()
        if not payload:
            self.statusBar().showMessage("Clipboard is empty — select blocks and press Ctrl+C first.", 4000)
            return
        if self.draw_action.isChecked():
            self.set_mode('select')
        self._paste_count += 1
        new_items = self.scene.paste_payload(payload, position, step=self._paste_count)
        self.statusBar().showMessage(f"Pasted {sum(isinstance(item, BlockItem) for item in new_items)} block(s) "
                                     f"and {sum(isinstance(item, LineItem) for item in new_items)} line(s).", 5000)

    def rotate_selected(self, clockwise=True):
        if self._text_editing():
            return
        selected_blocks = [item for item in self.scene.selectedItems() if isinstance(item, BlockItem)]
        if not selected_blocks:
            self.statusBar().showMessage("Select one or more blocks to rotate.", 4000)
            return
        for block in selected_blocks:
            self.scene.rotate_block(block, clockwise)

    def set_show_data_labels(self, visible):
        self.scene.set_show_labels(visible)
        self.scene.notify(physics_changed=False)

    # ---- undo / redo ----------------------------------------------------------------
    HISTORY_LIMIT = 200

    def _state_text(self):
        return json.dumps(self.scene.to_model(), sort_keys=True)

    def reset_history(self):
        """Start a fresh history at the current state (new / opened file)."""
        self._history_timer.stop()
        self._history = [self._state_text()]
        self._history_index = 0
        self._update_history_actions()

    def _request_history_record(self):
        if not self._restoring_history:
            self._history_timer.start()

    def _record_history(self):
        """Add the current state if it differs from the one at the history cursor
        (changes arriving together, e.g. a multi-block rotate, become one step)."""
        if self._restoring_history or not self._history:
            return
        state_text = self._state_text()
        if state_text == self._history[self._history_index]:
            return
        del self._history[self._history_index + 1:]
        self._history.append(state_text)
        if len(self._history) > self.HISTORY_LIMIT:
            del self._history[0]
        self._history_index = len(self._history) - 1
        self._update_history_actions()

    def _update_history_actions(self):
        self.undo_action.setEnabled(self._history_index > 0)
        self.redo_action.setEnabled(self._history_index < len(self._history) - 1)

    def undo(self):
        if self._text_editing():
            return
        if self._history_timer.isActive():          # record a change still waiting to be recorded
            self._history_timer.stop()
            self._record_history()
        if self._history_index > 0:
            self._restore_history(self._history_index - 1)
            self.statusBar().showMessage("Undone.", 2000)

    def redo(self):
        if self._text_editing():
            return
        if self._history_timer.isActive():
            self._history_timer.stop()
            self._record_history()
        if self._history_index < len(self._history) - 1:
            self._restore_history(self._history_index + 1)
            self.statusBar().showMessage("Redone.", 2000)

    def _restore_history(self, index):
        self._restoring_history = True
        try:
            self._history_index = index
            self.scene.load_model(json.loads(self._history[index]))
            self.labels_action.setChecked(self.scene.show_labels)
            self.scenario_panel.refresh()
            self.solve_now()
            self.set_modified(True)
        finally:
            self._history_timer.stop()
            self._restoring_history = False
        self._update_history_actions()

    def size_orifice(self, orifice_identifier):
        OrificeSizingDialog(self, orifice_identifier).exec()

    def find_lowest_speed(self, pump_identifier):
        LowestSpeedDialog(self, pump_identifier).exec()

    # ---- files ----------------------------------------------------------------------
    def _confirm_discard(self):
        if not self.modified:
            return True
        answer = QMessageBox.question(self, APPLICATION_NAME, "Discard unsaved changes?",
                                      QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        return answer == QMessageBox.StandardButton.Yes

    def new_file(self):
        if self._confirm_discard():
            self.load_model({'version': FILE_FORMAT_VERSION, 'blocks': [], 'lines': [], 'scenarios': []}, None)

    def open_dialog(self):
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open schematic", SCRIPT_DIRECTORY,
                                              "Pylon schematic (*.pylon.json *.prop.json *.json);;All files (*)")
        if path:
            self.open_file(path)

    def open_file(self, path):
        try:
            with open(path, encoding='utf-8') as file_handle:
                model = json.load(file_handle)
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Open", f"Could not open {path}:\n{error}")
            return
        self.load_model(model, path)

    def save(self):
        if not self.file_path:
            return self.save_as()
        try:
            with open(self.file_path, 'w', encoding='utf-8') as file_handle:
                json.dump(self.scene.to_model(), file_handle, indent=1)
        except OSError as error:
            QMessageBox.warning(self, "Save", str(error))
            return False
        self.set_modified(False)
        self.statusBar().showMessage(f"Saved {self.file_path}", 5000)
        return True

    def save_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save schematic",
                                              self.file_path or os.path.join(SCRIPT_DIRECTORY, 'schematic.pylon.json'),
                                              "Pylon schematic (*.pylon.json)")
        if not path:
            return False
        self.file_path = path
        return self.save()

    def export_image(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export image", os.path.join(SCRIPT_DIRECTORY, 'schematic.png'),
                                              "PNG image (*.png)")
        if not path:
            return
        self.scene.clearSelection()
        bounds = QRectF()
        for item in self.scene.items():
            if item.isVisible() and item is not self.scene.preview:
                bounds = bounds.united(item.sceneBoundingRect())
        bounds = bounds.adjusted(-30, -30, 30, 30)
        image = QImage(int(bounds.width() * 2), int(bounds.height() * 2), QImage.Format.Format_ARGB32)
        image.fill(QColor(THEME['background_deep']))
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.scene.render(painter, QRectF(image.rect()), bounds)
        painter.end()
        image.save(path)
        self.statusBar().showMessage(f"Exported {path}", 6000)

    def closeEvent(self, event):
        if self._confirm_discard():
            event.accept()
        else:
            event.ignore()


def main():
    application = QApplication(sys.argv)
    application.setApplicationName(APPLICATION_NAME)
    apply_application_theme(application)
    arguments = [argument for argument in sys.argv[1:] if not argument.startswith('-')]
    window = MainWindow(arguments[0] if arguments else None)
    window.show()
    sys.exit(application.exec())


if __name__ == '__main__':
    main()
