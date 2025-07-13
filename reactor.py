# Fission Reactor class
import math
from typing import Optional

from components import Component, FuelRod, ControlRod, CoolantChannel

R: float = 8.31446261815324
STANDARD_PRESSURE: float = 101325
ROOM_TEMPERATURE: float = 273
AIR_BOILING_POINT: float = 78.8
ZIRCALOY_HYDROGEN_REACTION_TEMPERATURE: float = 1500
THERMAL_CONDUCTIVITY: float = 45
WALL_THICKNESS: float = .1
COOLANT_WALL_THICKNESS: float = .06
SPECIFIC_HEAT_CAPACITY: float = 420
CONVECTIVE_HEAT_TRANSFER_COEFFICIENT: float = 10
POWER_DEFECT_COEFFICIENT: float = .016
DECAY_PRODUCT_RATE: float = .997
POISON_FRACTION: float = .063
CROSS_SECTION_RATIO: float = 4

class FissionReactor(object):
    reactor_layout: list[list[Optional[Component]]]
    fuel_rods: list[FuelRod]
    control_rods: list[ControlRod]
    coolant_channels: list[CoolantChannel]
    effective_control_rods: list[ControlRod]
    effective_coolant_channels: list[CoolantChannel]

    k: float
    control_rod_factor: float
    k_eff: float

    control_rod_insertion: float
    reactor_depth: int
    reactor_radius: float
    moderator_tipped: bool

    power: float # megawatts

    temperature: float = ROOM_TEMPERATURE
    pressure: float = STANDARD_PRESSURE
    exterior_pressure: float = STANDARD_PRESSURE

    coolant_boiling_point_standard_pressure: float
    coolant_exit_temperature: float
    prev_temperature: float

    coolant_heat_of_vaporization: float
    coolant_base_temperature: float
    fuel_depletion: float = -1
    neutron_poison_amount: float
    decay_products_amount: float
    env_temperature = ROOM_TEMPERATURE
    accumulated_hydrogen: float
    weighted_generation_time: float = 2
    max_temperature: float = 2000
    max_pressure: float = 15000000
    max_power: float = 3
    surface_area: float

    decay_neutrons: float
    neutron_flux: float
    neutron_to_power_conversion: float
    coolant_mass: float
    fuel_mass: float
    structural_mass: float
    control_rod_regulation_on: bool = True
    is_on: bool = False

    def __init__(self):
        pass

    def prepare_thermal_properties(self):
        pass

    def compute_geometry(self):
        pass

    def update_temperature(self):
        pass

    def update_pressure(self):
        pass

    def update_neutron_poisoning(self):
        pass

    def get_total_decay_neutrons(self) -> float:
        return self.neutron_poison_amount * 0.05 + self.decay_products_amount * 0.1 + self.decay_neutrons

    def update_power(self):
        pass

    def check_for_meltdown(self) -> bool:
        return self.temperature > self.max_temperature

    def check_for_explosion(self) -> bool:
        return self.pressure > self.max_pressure

    def add_component(self, component: Component, x: int, y: int) -> None:
        self.reactor_layout[x][y] = component

    def update_control_rod_insertion(self, insertion: float) -> None:
        self.control_rod_insertion = max(0.001, insertion)
        self.control_rod_factor = ControlRod.control_rod_factor(self.effective_control_rods, self.control_rod_insertion)

    def regulate_control_rods(self) -> None:
        if not self.is_on or not self.control_rod_regulation_on:
            return

        prev_insertion = self.control_rod_insertion
        if self.pressure > self.max_pressure * 0.8 \
                or self.temperature > (self.coolant_exit_temperature + self.max_temperature) / 2 \
                or self.temperature > self.max_temperature - 150 \
                or self.temperature - self.prev_temperature > 30:
            if self.k_eff > 0.99:
                self.control_rod_insertion += .004
        elif self.temperature > self.coolant_exit_temperature * 0.3 + self.coolant_base_temperature * 0.7:
            if self.k_eff > 1.01:
                self.control_rod_insertion += .008
            elif self.k_eff < 1.005:
                self.control_rod_insertion -= .001
        elif self.temperature > self.coolant_exit_temperature * 0.1 + self.coolant_base_temperature * 0.9:
            if self.k_eff > 1.025:
                self.control_rod_insertion += .012
            elif self.k_eff < 1.015:
                self.control_rod_insertion -= .004
        else:
            if self.k_eff > 1.1:
                self.control_rod_insertion += .02
            elif self.k_eff < 1.05:
                self.control_rod_insertion -= .006

        if prev_insertion != self.control_rod_insertion:
            self.control_rod_insertion = max(0.0, min(1.0, self.control_rod_insertion))
            self.control_rod_factor = ControlRod.control_rod_factor(self.effective_control_rods, self.control_rod_insertion)

    def turn_off(self):
        self.is_on = False
        self.max_power = 0
        self.k = 0
        self.k_eff = 0
        self.coolant_mass = 0
        self.fuel_mass = 0
        for i in range(len(self.reactor_layout)):
            for j in range(len(self.reactor_layout[i])):
                self.reactor_layout[i][j] = None
        self.fuel_rods.clear()
        self.control_rods.clear()
        self.coolant_channels.clear()
        self.effective_control_rods.clear()
        self.effective_coolant_channels.clear()
