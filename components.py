import math

class Fluid(object):
    temperature: float

class CoolantStats(object):
    hot_coolant: Fluid
    specific_heat_capacity: float
    moderator_factor: float
    cooling_factor: float
    boiling_point: float
    heat_of_vaporization: float
    accumulates_hydrogen: bool
    mass: float

COOLANTS = {}
COOLANTS_INVERSE = {}

def register_coolant(fluid: Fluid, coolant: CoolantStats):
    COOLANTS[fluid] = coolant
    COOLANTS_INVERSE[coolant] = fluid

def original_fluid(coolant: CoolantStats) -> Fluid:
    return COOLANTS_INVERSE[coolant]

class FuelStats(object):
    max_temperature: float
    duration: float
    slow_neutron_capture_cross_section: float
    fast_neutron_capture_cross_section: float
    slow_neutron_fission_cross_section: float
    fast_neutron_fission_cross_section: float
    neutron_generation_time: float
    released_neutrons: float
    required_neutrons: float = 1
    released_heat_energy: float
    decay_rate: float
    id: str

    def get_neutron_generation_time_category(self) -> int:
        if self.neutron_generation_time > 2:
            return 0
        elif self.neutron_generation_time > 1.25:
            return 1
        elif self.neutron_generation_time > 0.9:
            return 2
        else:
            return 3

    def get_fast_fission_multiplier(self) -> float:
        return self.fast_neutron_fission_cross_section * self.released_neutrons / self.required_neutrons

    def get_slow_fission_multiplier(self) -> float:
        return self.slow_neutron_fission_cross_section * self.released_neutrons / self.required_neutrons


class Component(object):
    moderation_factor: float
    max_temperature: float
    thermal_conductivity: float
    mass: float
    x: int
    y: int
    is_valid: bool
    index: int

    def __init__(self, moderation_factor: float, max_temperature: float, thermal_conductivity: float,
                 mass: float, is_valid: bool):
        self.moderation_factor = moderation_factor
        self.max_temperature = max_temperature
        self.thermal_conductivity = thermal_conductivity
        self.mass = mass
        self.is_valid = is_valid
        self.index = -1

    def set_pos(self, x: int, y: int) -> None:
        self.x = x
        self.y = y

    def get_absorption_factor(self, controls_inserted: bool, is_thermal: bool) -> float:
        return 0

    def same_position_as(self, component: "Component") -> bool:
        return component.x == self.x and component.y == self.y

    def get_distance(self, component: "Component") -> float:
        return math.sqrt((self.x - component.x) ** 2 + (self.y - component.y) ** 2)

class CoolantChannel(Component):
    coolant: CoolantStats
    weight: float
    related_fuel_rod_pairs: int
    partial_coolant: float

    def __init__(self, max_temperature: float, thermal_conductivity: float, coolant: CoolantStats, mass: float):
        super().__init__(coolant.moderator_factor, max_temperature, thermal_conductivity, mass, True)
        self.coolant = coolant
        self.weight = 0

    @staticmethod
    def normalize_weights(effective_channels: list["CoolantChannel"]) -> None:
        ch_sum = 0
        for channel in effective_channels:
            ch_sum += channel.weight
        for channel in effective_channels:
            channel.weight /= ch_sum

    def add_fuel_rod_pair(self) -> None:
        self.related_fuel_rod_pairs += 1

    def compute_weight_from_fuel_rod_map(self) -> None:
        self.weight = self.related_fuel_rod_pairs * 2


class FuelRod(Component):
    _fuel: FuelStats
    weight: float

    def __init__(self, max_temperature: float, thermal_conductivity: float, fuel: FuelStats, mass: float):
        super().__init__(0, max_temperature, thermal_conductivity, mass, True)
        self._fuel = fuel
        self.weight = 1

    def get_duration(self) -> float:
        return self._fuel.duration

    def get_fuel(self) -> FuelStats:
        return self._fuel

    def set_fuel(self, fuel: FuelStats) -> None:
        self._fuel = fuel
        self.max_temperature = fuel.max_temperature

    def get_neutron_generation_time(self):
        return self._fuel.neutron_generation_time


class ControlRod(Component):
    weight: float
    tip_moderation: bool
    related_fuel_rod_pairs: int

    def __init__(self, max_temperature: float, tip_moderation: bool, thermal_conductivity: float, mass: float):
        super().__init__(0, max_temperature, thermal_conductivity, mass, True)
        self.tip_moderation = tip_moderation
        self.weight = 0

    @staticmethod
    def normalize_weights(effective_control_rods: list["ControlRod"], total_weight: float, total_worth: float) -> None:
        for control in effective_control_rods:
            control.weight = control.weight / total_weight * total_worth

    @staticmethod
    def control_rod_factor(effective_control_rods: list["ControlRod"], insertion: float) -> float:
        crf = 0
        for rod in effective_control_rods:
            if rod.tip_moderation:
                if insertion <= 0.3:
                    crf -= insertion / 3 * rod.weight
                else:
                    crf -= (-11 / 7 * (insertion - 0.3) + 0.1) * rod.weight
            else:
                crf += insertion * rod.weight
        return crf

    def get_absorption_factor(self, controls_inserted: bool, is_thermal: bool) -> float:
        return 4 if controls_inserted else 0

    def add_fuel_rod_pair(self) -> None:
        self.related_fuel_rod_pairs += 1

    def compute_weight_from_fuel_rod_map(self) -> None:
        self.weight = self.related_fuel_rod_pairs * 4
