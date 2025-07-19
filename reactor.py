# Fission Reactor class
import math
from os import MFD_ALLOW_SEALING

import vecmath
from vecmath import *
from typing import Optional, Callable

from components import Component, FuelRod, ControlRod, CoolantChannel, CoolantStats, original_fluid

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
        id_rod = 0
        id_control = 0
        id_channel = 0
        for i in range(len(self.reactor_layout)):
            for j in range(len(self.reactor_layout[i])):
                comp = self.reactor_layout[i][j]
                comp.set_pos(i, j)
                self.max_temperature = min(self.max_temperature, comp.max_temperature)
                self.structural_mass += comp.mass
                if isinstance(comp, FuelRod):
                    comp.index = id_rod
                    self.fuel_rods.append(comp)
                    id_rod += 1
                elif isinstance(comp, ControlRod):
                    comp.index = id_control
                    self.control_rods.append(comp)
                    id_control += 1
                elif isinstance(comp, CoolantChannel):
                    comp.index = id_channel
                    self.coolant_channels.append(comp)
                    id_channel += 1

    def compute_k(self, add_to_effective_lists: bool, control_rods_inserted: bool):
        geometric_matrix_neutrons = [[0.0] * len(self.fuel_rods)] * len(self.fuel_rods)

        for i in range(len(self.fuel_rods)):
            for j in range(i):
                mij = 0.0
                saij = 0.0
                faij = 0.0
                rod_one = self.fuel_rods[i]
                rod_two = self.fuel_rods[j]

                prev_x = self.fuel_rods[i].x
                prev_y = self.fuel_rods[i].y
                resolution = 100
                for t in range(resolution):
                    x = round((rod_two.x - rod_one.x) * (t / resolution) + self.fuel_rods[i].x)
                    y = round((rod_two.y - rod_one.y) * (t / resolution) + self.fuel_rods[i].y)
                    if x < 0 or x > len(self.reactor_layout) - 1 or y < 0 or y > len(self.reactor_layout) - 1:
                        continue

                    component = self.reactor_layout[x][y]
                    if component is None:
                        continue

                    if component.moderation_factor > 0:
                        mij += component.moderation_factor
                        saij = (faij + saij) / 2

                    if (not component.same_position_as(self.fuel_rods[i])
                            and not component.same_position_as(self.fuel_rods[j])):
                        saij += component.get_absorption_factor(control_rods_inserted, True)
                        faij += component.get_absorption_factor(control_rods_inserted, False)

                    if not add_to_effective_lists or (x == prev_x and y == prev_y):
                        continue

                    prev_x = x
                    prev_y = y

                    if isinstance(component, CoolantChannel) or isinstance(component, ControlRod):
                        component.add_fuel_rod_pair()

                mij /= resolution
                faij /= resolution
                saij /= resolution

                dist = rod_one.get_distance(rod_two)
                unabsorbed_fast = math.exp(-faij * dist) / dist
                unabsorbed_slow = math.exp(-saij * dist) / dist
                fast = math.exp(-mij * dist) / dist
                slow = (1 / dist - fast) * unabsorbed_slow
                fast *= unabsorbed_fast

                slow_neutron_fission_mult = rod_two.get_fuel().get_slow_fission_multiplier()
                fast_neutron_fission_mult = rod_two.get_fuel().get_fast_fission_multiplier()
                geometric_matrix_neutrons[i][j] = slow * slow_neutron_fission_mult + fast * fast_neutron_fission_mult

                slow_neutron_fission_mult = rod_one.get_fuel().get_slow_fission_multiplier()
                fast_neutron_fission_mult = rod_one.get_fuel().get_fast_fission_multiplier()
                geometric_matrix_neutrons[j][i] = slow * slow_neutron_fission_mult + fast * fast_neutron_fission_mult

        vector = [1.0] * len(self.fuel_rods)
        for i in range(10):
            vector = vecmath.normalize(vector)
            vector = vecmath.multiply(geometric_matrix_neutrons, vector)

        k_calc = vecmath.magnitude(vector)
        if add_to_effective_lists:
            vector = vecmath.linear_normalize(vector)
            for i in range(len(self.fuel_rods)):
                self.fuel_rods[i].weight = vector[i]

        k_calc *= self.reactor_depth / (1 + self.reactor_depth)
        return k_calc

    def compute_geometry(self):
        self.effective_control_rods.clear()
        self.effective_coolant_channels.clear()
        self.moderator_tipped = False

        self.k = self.compute_k(True, False)
        k_experimental = self.compute_k(False, True)

        self.compute_control_rod_weights(((self.k - 1) / self.k) - ((k_experimental - 1) / k_experimental))

        self.neutron_to_power_conversion = 0
        self.decay_neutrons = 0

        for i_idx, i in enumerate(self.fuel_rods):
            self.neutron_to_power_conversion += i.get_fuel().released_heat_energy / i.get_fuel().released_neutrons
            self.decay_neutrons += i.get_fuel().decay_rate

        if len(self.fuel_rods) > 1:
            self.neutron_to_power_conversion /= len(self.fuel_rods)
            self.max_power = self.calculate_max_power()
        else:
            self.k = 0.00001
            self.max_power = 0.1 * 0.1

        self.compute_coolant_channel_weights()
        self.control_rod_factor = ControlRod.control_rod_factor(self.effective_control_rods, self.control_rod_insertion)
        self.prepare_initial_conditions()

    def compute_control_rod_weights(self, total_worth: float):
        total_weight = 0
        for rod in self.control_rods:
            rod.compute_weight_from_fuel_rod_map()
            if rod.weight > 0:
                self.effective_control_rods.append(rod)
                total_weight += rod.weight

        ControlRod.normalize_weights(self.effective_control_rods, total_weight, total_worth)

    def compute_coolant_channel_weights(self):
        for channel in self.coolant_channels:
            channel.compute_weight_from_fuel_rod_map()
            if channel.weight > 0:
                self.effective_coolant_channels.append(channel)

        CoolantChannel.normalize_weights(self.effective_coolant_channels)

    def reset_fuel_depletion(self):
        self.fuel_depletion = 0

    def prepare_initial_conditions(self):
        self.coolant_base_temperature = 0
        self.coolant_boiling_point_standard_pressure = 0
        self.coolant_exit_temperature = 0
        self.coolant_heat_of_vaporization = 0
        self.weighted_generation_time = 0

        for rod in self.fuel_rods:
            self.weighted_generation_time += rod.get_neutron_generation_time()
        self.weighted_generation_time /= len(self.fuel_rods)

        for channel in self.coolant_channels:
            prop = channel.coolant
            fluid = original_fluid(prop)
            if fluid is not None:
                self.coolant_base_temperature += fluid.temperature
            self.coolant_boiling_point_standard_pressure += prop.boiling_point
            self.coolant_exit_temperature += prop.hot_coolant.temperature
            self.coolant_heat_of_vaporization += prop.heat_of_vaporization

        if len(self.coolant_channels) > 0:
            self.coolant_base_temperature /= len(self.coolant_channels)
            self.coolant_boiling_point_standard_pressure /= len(self.coolant_channels)
            self.coolant_exit_temperature /= len(self.coolant_channels)
            self.coolant_heat_of_vaporization /= len(self.coolant_channels)

            if self.coolant_base_temperature == 0:
                self.coolant_base_temperature = self.env_temperature
            if self.coolant_boiling_point_standard_pressure == 0:
                self.coolant_boiling_point_standard_pressure = AIR_BOILING_POINT

        self.is_on = True

    def make_coolant_flow(self) -> float:
        heat_removed = 0
        self.coolant_mass = 0

        for channel in self.coolant_channels:
            # this tries to drain up to 16,000 and then "drained" is the actual amount we could drain
            # if the input hatch has less, then it would be less.
            drained = 16000
            prop = channel.coolant
            coolant_temp = original_fluid(prop).temperature

            cooled_temp = prop.hot_coolant.temperature
            if cooled_temp > self.temperature:
                continue

            heat_removed_per_liter = prop.specific_heat_capacity / 14 * (cooled_temp - coolant_temp)

            heat_flux_per_area_and_temp = 1 / (1 / prop.cooling_factor + COOLANT_WALL_THICKNESS / THERMAL_CONDUCTIVITY)
            ideal_heat_flux = heat_flux_per_area_and_temp * 4 * self.reactor_depth * (self.temperature - cooled_temp)

            ideal_fluid_used = ideal_heat_flux / heat_removed_per_liter
            capped_fluid_used = min(drained, int(ideal_fluid_used))

            remaining_space = 16000 # amount of space in the output hatch
            actual_flow_rate = min(remaining_space, int(capped_fluid_used + channel.partial_coolant))
            channel.partial_coolant += capped_fluid_used - actual_flow_rate

            # actually drain actual_flow_rate from the input handler, and add an equivalent amount of
            # hot coolant to the output handler

            if prop.accumulates_hydrogen and self.temperature > ZIRCALOY_HYDROGEN_REACTION_TEMPERATURE:
                boiling_point = self.coolant_boiling_point_standard_pressure
                if self.temperature > boiling_point:
                    self.accumulated_hydrogen += (self.temperature - boiling_point) / boiling_point
                elif actual_flow_rate < min(remaining_space, int(ideal_fluid_used)):
                    self.accumulated_hydrogen += ((self.temperature - ZIRCALOY_HYDROGEN_REACTION_TEMPERATURE)
                                                  / ZIRCALOY_HYDROGEN_REACTION_TEMPERATURE)

            self.coolant_mass += capped_fluid_used * prop.mass
            heat_removed += capped_fluid_used * heat_removed_per_liter

        self.coolant_mass /= 1000
        self.accumulated_hydrogen *= .98
        return heat_removed

    def calculate_max_power(self) -> float:
        hypothetical_temp = min(self.max_temperature, ZIRCALOY_HYDROGEN_REACTION_TEMPERATURE)
        heat_removed = 0
        for channel in self.coolant_channels:
            prop = channel.coolant
            coolant_temp = original_fluid(prop).temperature

            cooled_temp = prop.hot_coolant.temperature
            if cooled_temp > hypothetical_temp:
                continue

            heat_removed_per_liter = prop.specific_heat_capacity / 14 * (cooled_temp - coolant_temp)

            heat_flux_per_area_and_temp = 1 / (1 / prop.cooling_factor + COOLANT_WALL_THICKNESS / THERMAL_CONDUCTIVITY)
            ideal_heat_flux = heat_flux_per_area_and_temp * 4 * self.reactor_depth * (hypothetical_temp - cooled_temp)

            ideal_fluid_used = ideal_heat_flux / heat_removed_per_liter
            heat_removed += ideal_fluid_used * heat_removed_per_liter

        time_constant = SPECIFIC_HEAT_CAPACITY * (1 / CONVECTIVE_HEAT_TRANSFER_COEFFICIENT +
                                                  WALL_THICKNESS / THERMAL_CONDUCTIVITY) / self.surface_area

        return ((hypothetical_temp - self.env_temperature) *
                (time_constant * (self.coolant_mass + self.structural_mass + self.fuel_mass)) + heat_removed) / 1e6

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
