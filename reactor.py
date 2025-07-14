# Fission Reactor class
import math
from vecmath import *
from typing import Optional, Callable

from components import Component, FuelRod, ControlRod, CoolantChannel, CoolantStats

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

def response_function(target: float, current: float, critical_rate: float) -> float:
    if current < 0:
        if critical_rate < 1:
            return 0
        else:
            current = 0.1
    exp_decay = math.exp(-critical_rate)
    return current * exp_decay + target * (1 - exp_decay)

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

    custom_regulation_logic: Optional[Callable[["FissionReactor"], float]] = None

    def __init__(self, size: int, depth: int, insertion: float):
        self.reactor_layout = [[None] * size] * size
        self.reactor_depth = depth
        self.reactor_radius = float(size) / 2 + 1.5
        self.control_rod_insertion = max(0.001, insertion)
        self.fuel_rods = []
        self.control_rods = []
        self.coolant_channels = []
        self.effective_control_rods = []
        self.effective_coolant_channels = []
        self.surface_area = (self.reactor_radius ** 2 * math.pi * 2 +
                             self.reactor_depth * self.reactor_radius * math.pi * 2)
        self.structural_mass = self.reactor_depth * self.reactor_radius ** 2 * math.pi * 300

    def response_function_temperature(self, env_temp: float, current_temp: float, heat_added: float,
                                      heat_absorbed: float) -> float:
        current_temp = max(0.1, current_temp)
        heat_absorbed = max(0.0, heat_absorbed)
        time_constant = (SPECIFIC_HEAT_CAPACITY *
                        (1 / CONVECTIVE_HEAT_TRANSFER_COEFFICIENT + WALL_THICKNESS / THERMAL_CONDUCTIVITY) /
                        self.surface_area)
        exp_decay = math.exp(-time_constant)
        effective_env_temp = (env_temp + (heat_added - heat_absorbed) /
                             (time_constant * (self.coolant_mass + self.structural_mass + self.fuel_mass)))
        return current_temp * exp_decay + effective_env_temp * (1 - exp_decay)

    def prepare_thermal_properties(self):
        pass

    def compute_geometry(self):
        pass

    def make_coolant_flow(self) -> float:
        pass

    def calculate_max_power(self):
        pass

    def coolant_boiling_point(self, coolant: Optional[CoolantStats] = None) -> float:
        if coolant is None:
            return self.coolant_boiling_point_standard_pressure
        else:
            bp = coolant.boiling_point
            return bp if bp != 0 else self.coolant_boiling_point_standard_pressure

    def update_temperature(self):
        self.prev_temperature = self.temperature
        self.temperature = self.response_function_temperature(self.env_temperature, self.temperature,
                                                              self.power * 1e6, 0)
        self.temperature = min(self.max_temperature, self.temperature)
        heat_removed = self.make_coolant_flow()
        self.temperature = self.response_function_temperature(self.env_temperature, self.prev_temperature,
                                                              self.power * 1e6, heat_removed)
        self.temperature = max(self.temperature, self.coolant_base_temperature)

    def update_pressure(self):
        use_exterior = self.temperature > self.coolant_boiling_point() and self.is_on
        self.pressure = response_function(self.exterior_pressure if use_exterior else 1000. * R * self.temperature,
                                          self.pressure, 0.2)

    def update_neutron_poisoning(self):
        self.neutron_poison_amount += self.decay_products_amount * (1 - DECAY_PRODUCT_RATE) * POISON_FRACTION
        self.neutron_poison_amount *= DECAY_PRODUCT_RATE * math.exp(-CROSS_SECTION_RATIO * self.power / self.surface_area)

    def get_total_decay_neutrons(self) -> float:
        return self.neutron_poison_amount * 0.05 + self.decay_products_amount * 0.1 + self.decay_neutrons

    def update_power(self):
        if self.is_on:
            self.neutron_flux += self.get_total_decay_neutrons()
            self.k_eff = 1 / ((1 / self.k) + POWER_DEFECT_COEFFICIENT * (self.power / self.max_power) +
                              self.neutron_poison_amount * CROSS_SECTION_RATIO /
                              self.surface_area + self.control_rod_factor)
            self.k_eff = max(0., self.k_eff)

            inverse_reactor_period = (self.k_eff - 1) / self.weighted_generation_time
            self.neutron_flux *= math.exp(inverse_reactor_period)

            self.fuel_depletion += self.neutron_flux * self.reactor_depth
            self.decay_products_amount += max(self.neutron_flux, 0.) / 1000

            self.power = self.neutron_flux * self.neutron_to_power_conversion
        else:
            self.neutron_flux *= 0.5
            self.power *= 0.5
        self.decay_products_amount *= DECAY_PRODUCT_RATE

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
        if self.custom_regulation_logic is not None:
            self.control_rod_insertion = self.custom_regulation_logic(self)
        elif self.pressure > self.max_pressure * 0.8 \
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
