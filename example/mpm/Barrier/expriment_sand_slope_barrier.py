# ls mkdir pwd cd ps
# vim: exit, edit file, save file
# git: set up a repository on github; git add, git commit -m, git push, git pull, git log
import sys
import numpy as np
import taichi as ti

sys.path.append('c:/Users/86198/Desktop/Geo/GeoTaichi')
from geotaichi import *

init()

mpm = MPM()

# Flume and simulation parameters
theta = np.radians(45)  # Flume inclination angle
gravity_x = 9.8 * np.sin(theta)
gravity_z = -9.8 * np.cos(theta)
flume_length = 2.4
flume_width = 0.3
flume_height = 0.5# Height sufficient to cover sand deposit and wall

# Set the simulation configuration
mpm.set_configuration(domain=ti.Vector([flume_length, flume_width, flume_height]),
                      background_damping=0.0,
                      gravity=ti.Vector([-6.929646456, 0., -6.929646456]),
                      alphaPIC=0.00,
                      mapping="Newmark",
                      shape_function="GIMP",
                      stabilize=None)

# Solver configuration
mpm.set_solver(solver={
    "Timestep": 1e-5,
    "SimulationTime": 0.50,
    "SaveInterval": 0.01
})

# Memory allocation
mpm.memory_allocate(memory={
    "max_material_number": 2,
    "max_particle_number": 1e6,
    "max_constraint_number": {
        "max_velocity_constraint": 5e7,
         "max_friction_constraint":   5e4
    }
})

# Material properties
density_sand = 1379
youngs_modulus_sand = 21.6e6
poissons_ratio_sand = 0.3
density_wall = 2300
youngs_modulus_wall = 40e9
poissons_ratio_wall = 0.15

# Sand material
mpm.add_material(model="DruckerPrager",
                 material={
                     "MaterialID": 1,
                     "Density": density_sand,
                     "YoungModulus": youngs_modulus_sand,
                     "PoissonRatio": poissons_ratio_sand,
                     "Friction": np.deg2rad(30),
                     "Dilation": 0,
                     "Cohesion": 0
                 })


# Sand deposit region
sand_deposit_height = 0.3
mpm.add_region(region={
    "Name": "sand_deposit",
    "Type": "Rectangle",
    "BoundingBoxPoint": ti.Vector([1.9, 0., 0.01]),
    "BoundingBoxSize": ti.Vector([0.5, 0.3, 0.3]),
    "zdirection": ti.Vector([0., 0., 1.])
})

mpm.add_region(region={
    "Name": "Slope",
    "Type": "Rectangle",
    # 确保其他属性正确定义
    "BoundingBoxPoint": ti.Vector([0, 0., 0]),
    "BoundingBoxSize": ti.Vector([2.4, 0.3, 0.01]),
    "zdirection": ti.Vector([0, 0, 1])
})

# Element configuration
element_size = 0.01
mpm.add_element(element={
    "ElementType": "R8N3D",
    "ElementSize": ti.Vector([element_size, element_size, element_size]),
    "Contact": {
        "ContactDetection": True,
        "Friction": 0.3,
        "CutOff": 0.8
    }
})

# Add bodies for sand
mpm.add_body(body={
    "Template": [
        {
            "RegionName": "sand_deposit",
            "nParticlesPerCell": 2,
            "BodyID": 0,
            "MaterialID": 1,
            "ParticleStress": {
                "InternalStress": ti.Vector([0., 0., 0., 0., 0., 0.])
            },
            "InitialVelocity": ti.Vector([0., 0., 0.]),
            "FixVelocity": ["Free", "Free", "Free"]
        },
        {
            "RegionName": "Slope",
            "nParticlesPerCell": 1,
            "BodyID": 1,
            "RigidBody": True,
            "ParticleStress": {
                "GravityField": False,
                "InternalStress": ti.Vector([0., 0., 0., 0., 0., 0.])
            },
            "InitialVelocity": ti.Vector([0., 0., 0.]),
            "FixVelocity": ["Fix", "Fix", "Fix"]
        }
    ]
})


# Example of boundary conditions setup
# 设定边界条件列表
boundary_conditions = [
    # 斜坡的摩擦约束
    {
        "BoundaryType": "FrictionConstraint",
        "Friction": 0.3,  # 设定摩擦系数，根据需要调整
        "StartPoint": [0., 0., 0.01],  # 根据你的斜坡region的BoundingBoxPoint调整
        "EndPoint": [2.4, 0.3, 0],  # 根据你的斜坡region的BoundingBoxSize调整
        "Norm": [0, 0, 1],  # 假设斜坡平行于z方向
    },
    # 墙的速度约束
    {
        "BoundaryType": "VelocityConstraint",
        "Velocity": [0, None, None],  # 墙是固定的，不动
        "StartPoint": [0., 0., 0.01],  # 根据你的墙region的BoundingBoxPoint调整
        "EndPoint": [0.1, 0.3, 0.3],  # 根据你的墙region的BoundingBoxSize调整
    },
    #Displacement constraint for 斜坡
    #{
        #"BoundaryType": "DisplacementConstraint",
        #"StartPoint": [0., 0., 0.01],  # 根据你的斜坡region的BoundingBoxPoint调整
        #"EndPoint": [2.4*0.1, 0.3*0.1, 0],  # 根据你的斜坡region的BoundingBoxSize调整
        #"Norm":[0,0,1],
    #}
    
]# 将边界条件列表添加到MPM仿真中
mpm.add_boundary_condition(boundary=boundary_conditions)


# Save data configuration
mpm.select_save_data(particle=True, grid=True)
# Run the simulation
mpm.run()
# Post-processing
mpm.postprocessing(start_file=0, end_file=50, write_background_grid=True)
