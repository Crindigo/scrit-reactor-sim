# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import reactor
from components import Fluid, CoolantStats, FuelStats, ControlRod, FuelRod, CoolantChannel

water_coolant = CoolantStats()
water_coolant.fluid = Fluid(293, "distilled_water")
water_coolant.hot_coolant = Fluid(558, "high_pressure_steam")
water_coolant.moderator_factor = 2
water_coolant.cooling_factor = 1000
water_coolant.boiling_point = 373
water_coolant.heat_of_vaporization = 2260000
water_coolant.specific_heat_capacity = 4168
water_coolant.accumulates_hydrogen = True
water_coolant.slow_absorption_factor = 0.1875
water_coolant.fast_absorption_factor = 0.0625
water_coolant.mass = 10 # 1*2 + 8

leu235 = FuelStats()
leu235.id = "LEU-235"
leu235.max_temperature = 1500
leu235.duration = 75000
leu235.neutron_generation_time = 3.5
leu235.slow_neutron_capture_cross_section = 1.8
leu235.slow_neutron_fission_cross_section = 1.8
leu235.required_neutrons = 1
leu235.released_neutrons = 2.5
leu235.released_heat_energy = 0.01
leu235.decay_rate = 0.025

def make_fuel_rod():
    return FuelRod(leu235.max_temperature, 1, leu235, 650)

def make_control_rod():
    return ControlRod(100000, False, 1, 800)

def make_coolant_channel():
    return CoolantChannel(100050, 0, water_coolant, 1000)

def run_reactor():
    # size=3 means max 3x3 of components
    # it is calculated as multiblock diameter - 2 (to account for walls)
    # so, valid values should be 3,5,7,9,11,13, even though on bigger sizes you can't use (0,0) for example.
    fr = reactor.FissionReactor(3, 1, 0.5)
    # fuel rods on 4 corners
    fr.add_component(make_fuel_rod(), 0, 0)
    fr.add_component(make_fuel_rod(), 2, 0)
    fr.add_component(make_fuel_rod(), 0, 2)
    fr.add_component(make_fuel_rod(), 2, 2)

    # control rods on vertical between fuel rods
    fr.add_component(make_control_rod(), 0, 1)
    fr.add_component(make_control_rod(), 2, 1)

    # line of coolant channels
    fr.add_component(make_coolant_channel(), 1, 0)
    fr.add_component(make_coolant_channel(), 1, 1)
    fr.add_component(make_coolant_channel(), 1, 2)

    fr.prepare_thermal_properties()
    fr.compute_geometry()

    print(repr(fr.__dict__))

    for i in range(1200):
        fr.update_power()
        fr.update_temperature()
        fr.update_pressure()
        fr.update_neutron_poisoning()
        fr.regulate_control_rods()
        if fr.check_for_meltdown():
            print("Meltdown")
            break

        if fr.check_for_explosion():
            print("Explosion")
            break

        if i % 10 == 0:
            print("%d | Keff = %.5f, Ctrl = %.2f, Pwr: %.2f, Temp: %.2f, Pr: %.2f"
                  % (i, fr.k_eff, fr.control_rod_insertion, fr.power, fr.temperature, fr.pressure))

def print_reactor(fr: reactor.FissionReactor):
    print("K = %.2f")

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_reactor()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
