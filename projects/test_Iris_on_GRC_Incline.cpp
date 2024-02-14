// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Ruochun Zhang, Luning Bakke
// =============================================================================
//
// Demo to show Rover operating on incline of GRC-1 simulant, with DEM-Engine
// providing the DEM simulation support
//
// =============================================================================

#include "chrono_models/robot/iris/Iris.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"

#include <DEM/API.h>
// #include <core/utils/DEMEPaths.hpp>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;
using namespace std::filesystem;

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::iris;

const double math_PI = 3.1415927;

IrisWheelType wheel_type = IrisWheelType::RealWheel;

// ChVector to/from float3
inline float3 ChVec2Float(const ChVector<>& vec) {
    return make_float3(vec.x(), vec.y(), vec.z());
}
inline ChVector<> Float2ChVec(float3 f3) {
    return ChVector<>(f3.x, f3.y, f3.z);
}
inline float4 ChQ2Float(const ChQuaternion<>& Q) {
    float4 f4;
    f4.w = Q.e0();
    f4.x = Q.e1();
    f4.y = Q.e2();
    f4.z = Q.e3();
    return f4;
}

// Load clump type, set clump material properties and load particle checkpoint file
void PrepareParticles(DEMSolver& DEMSim, std::shared_ptr<DEMMaterial> wheel_material);

void SaveParaViewFiles(Iris& rover, path& rover_dir, unsigned int frame_number);


// E, nu, CoR, mu, COR for soil parameters
float mu = 0.4;  // friction among soil particles
float mu_wheel = 0.8;  // friction among wheel and terrain
float mu_wall = 1.;
float CoR = 0.25;
float E = 1e8;


// Define the simulation world domain (small size for debug)
double world_x_size = 1.0;
double world_y_size = 1.0;
//    double world_x_size = 4.0f;
//    double world_y_size = 2.0f;


int main(int argc, char* argv[]) {
    SetChronoDataPath(CHRONO_DATA_DIR);
    SetDEMEDataPath(DEME_DATA_DIR);
    std::cout << "DEME dir is " << GetDEMEDataPath() << std::endl;

    // `World'
    float G_mag = 9.81;
    float step_size = 2e-6;  // 2e-6; // 1e-6 for 15 deg and above, perhaps

    // Use wheel geometry from the fact sheet
    float wheel_rad = 0.098;
    float wheel_width = 0.034;
    float wheel_mass = 0.028;
    float chassis_mass = 2.186;
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
    float3 wheel_MOI = make_float3(wheel_IXX, wheel_IYY, wheel_IXX);

    std::string wheel_obj_path = GetChronoDataFile("robot/iris/obj/iris_wheel.obj");

    std::cout << "wheel obj directory: " << wheel_obj_path << std::endl;
    // Create a Chrono::Engine physical system
    float Slope_deg;

    if (argc == 1)
            Slope_deg = 0;
    if (argc == 2) {
            Slope_deg = atof(argv[1]);
    }

    std::cout << "SLOPE IS "  << Slope_deg << std::endl;

    double G_ang = Slope_deg * math_PI / 180.;


    // Chrono system for iris wheel
    ChSystemNSC sys;
    // consistent gravity in DEM and Chrono
    ChVector<double> G = ChVector<double>(-G_mag * std::sin(G_ang), 0, -G_mag * std::cos(G_ang));
    sys.Set_G_acc(G);

    const int nW = 4;  // 4 wheels

    // Create the rover
    float w_r = 0.2;  // rover wheel driving speed (rad/sec)
    auto driver = chrono_types::make_shared<IrisSpeedDriver>(0.0, w_r);
    Iris iris(&sys, wheel_type);
    iris.SetSpeedDriver(driver);

    float z_offset = 0.08;
    // initial position of the rover
    iris.Initialize(ChFrame<>(ChVector<>(-0.2, -0.0, -0.225), QUNIT));

    // Get wheels and bodies to set up SCM patches
    std::vector<std::shared_ptr<ChBodyAuxRef>> Wheels;
    std::vector<ChVector<>> wheel_pos;
    Wheels.push_back(iris.GetWheel(IrisWheelID::LF)->GetBody());
    Wheels.push_back(iris.GetWheel(IrisWheelID::RF)->GetBody());
    Wheels.push_back(iris.GetWheel(IrisWheelID::LB)->GetBody());
    Wheels.push_back(iris.GetWheel(IrisWheelID::RB)->GetBody());

    auto Body_1 = iris.GetChassis()->GetBody();
    std::cout << "Rover mass: " << iris.GetRoverMass() << std::endl;
    std::cout << "Wheel mass: " << iris.GetWheelMass() << std::endl;

    for (int i = 0; i < nW; i++) {
        wheel_pos.push_back(Wheels[i]->GetFrame_REF_to_abs().GetPos());
    }

    //////////////////////////////////////////////
    // DEM setup
    //////////////////////////////////////////////

    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::XYZ);

    // Family 1 is for fixed ground which does not participate the force calc.
    DEMSim.SetFamilyFixed(1);
    DEMSim.DisableContactBetweenFamilies(1, 1);
    DEMSim.DisableContactBetweenFamilies(1, 255);

    DEMSim.SetCollectAccRightAfterForceCalc(true);
    auto mat_type_wheel = DEMSim.LoadMaterial({{"E", E}, {"nu", 0.3}, {"CoR", CoR}, {"mu", mu_wheel}, {"Crr", 0.00}});

    PrepareParticles(DEMSim, mat_type_wheel);

    ///////////////////////
    // Add wheel in DEM
    ///////////////////////

    // Instantiate this wheel
    std::cout << "Making wheels..." << std::endl;
    DEMSim.SetFamilyPrescribedAngVel(100);
    DEMSim.SetFamilyPrescribedLinVel(100);
    std::vector<std::shared_ptr<DEMTracker>> trackers;
    // std::vector<std::shared_ptr<DEMClumpBatch>> DEM_Wheels;
    std::vector<std::shared_ptr<DEMMeshConnected>> DEM_Wheels;
    for (int i = 0; i < nW; i++) {

        auto wheel_mesh = DEMSim.AddWavefrontMeshObject(wheel_obj_path, mat_type_wheel);
        DEM_Wheels.push_back(wheel_mesh);


        // If this is one of the left wheels
        if (i == 1 || i == 3) {
             DEM_Wheels[i]->Mirror(make_float3(0,0,0), make_float3(0,1,0));
         }
        DEM_Wheels[i]->SetFamily(100);
        DEM_Wheels[i]->SetMass(wheel_mass);
        DEM_Wheels[i]->SetMOI(wheel_MOI);
        trackers.push_back(DEMSim.Track(DEM_Wheels[i]));
    }
    DEMSim.DisableFamilyOutput(100);  // no need outputting wheels (if it's mesh, actually won't be outputted anyway)
    // DEMSim.DisableContactBetweenFamilies(100, 0);  // No wheel-ground contact while settling
    std::cout << "Total num of triangles in a wheel: " << DEM_Wheels[0]->GetNumTriangles() << std::endl;

    //////
    // Make ready for DEM simulation
    ///////
    // Some inspectors
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
    auto partial_mass_finder = DEMSim.CreateInspector("clump_mass", "return (Z <= -0.41);");
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(ChVec2Float(G));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEMSim.SetCDUpdateFreq(30);
    // DEMSim.SetInitBinSize(scales.at(1));
    DEMSim.SetInitBinNumTarget(5e6);

    DEMSim.SetExpandSafetyAdder(0.2);
    DEMSim.SetExpandSafetyMultiplier(1.);
    DEMSim.SetErrorOutVelocity(1e8);

    DEMSim.Initialize();
    for (const auto& tracker : trackers) {
        std::cout << "A tracker is tracking owner " << tracker->obj->ownerID << std::endl;
    }
    std::cout << "Time step size is " << step_size << std::endl;
    std::cout << "End initialization" << std::endl;

    float time_end = 15.0;
    unsigned int fps = 20;
    // unsigned int move_box_ps = 1;
    unsigned int report_freq = 2000;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    float frame_accu_thres = 1.0 / fps;
    unsigned int report_steps = (unsigned int)(1.0 / (report_freq * step_size));

    // path out_dir = current_path();
    path out_dir = GetChronoOutputPath();
    create_directory(out_dir);
    out_dir += "/Iris_thicker_wheel" + std::to_string(Slope_deg);
    path rover_dir = out_dir / "./rover";
    create_directory(out_dir);
    create_directory(rover_dir);
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    ///////////////////////////////////////////
    // Real simulation
    ///////////////////////////////////////////

    // Timers
    std::chrono::high_resolution_clock::time_point h_start, d_start;
    std::chrono::high_resolution_clock::time_point h_end, d_end;
    std::chrono::duration<double> h_total, d_total;
    h_total = std::chrono::duration<double>(0);
    d_total = std::chrono::duration<double>(0);

    std::vector<ChQuaternion<>> wheel_rot(4);
    std::vector<ChVector<>> wheel_vel(4);
    std::vector<ChVector<>> wheel_angVel(4);
    float max_v;
    int change_step = 0;
    float frame_accu = frame_accu_thres;

    // Put the wheels somewhere that won't affect simulation
    {
        trackers[0]->SetPos(make_float3(1, 0.5, 0.75));
        trackers[1]->SetPos(make_float3(1, -0.5, 0.75));
        trackers[2]->SetPos(make_float3(-1, 0.5, 0.75));
        trackers[3]->SetPos(make_float3(-1, -0.5, 0.75));
    }


    SaveParaViewFiles(iris, rover_dir, 0);

    // Settle first, then put the wheel in place, then let the wheel sink in initially
    for (float t = 0; t < 0.02; t += frame_accu_thres) {
        std::cout << "Num contacts: " << DEMSim.GetNumContacts() << std::endl;
        DEMSim.ShowThreadCollaborationStats();
        DEMSim.ShowTimingStats();
        DEMSim.DoDynamicsThenSync(frame_accu_thres);
    }

    unsigned int chrono_update_freq = 10;
    // Active box domain
    float box_halfsize_x = wheel_rad * 1.2;
    float box_halfsize_y = wheel_width * 2.;


    for (float t = 0; t < time_end; t += step_size, curr_step++, frame_accu += step_size) {


        if (curr_step % chrono_update_freq == 0) {
            for (int i = 0; i < nW; i++) {
                wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
                trackers[i]->SetPos(ChVec2Float(wheel_pos[i]));
                wheel_rot[i] = Wheels[i]->GetFrame_REF_to_abs().GetRot();
                trackers[i]->SetOriQ(ChQ2Float(wheel_rot[i]));
                wheel_vel[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos_dt();
                trackers[i]->SetVel(ChVec2Float(wheel_vel[i]));
                wheel_angVel[i] = Wheels[i]->GetFrame_REF_to_abs().GetWvel_par();
                trackers[i]->SetAngVel(ChVec2Float(wheel_angVel[i]));
            }
        }

        // at every frame step, write output and move active box domain
        if (frame_accu >= frame_accu_thres) {
            frame_accu = 0.;
            std::cout << "Frame: " << currframe << std::endl;
            std::cout << "Num contacts: " << DEMSim.GetNumContacts() << std::endl;
            std::cout << h_total.count() << " seconds spent on host" << std::endl;
            std::cout << d_total.count() << " seconds spent on device" << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            SaveParaViewFiles(iris, rover_dir, currframe);
            currframe++;
            DEMSim.ShowTimingStats();


            if (t > 0.2) {
                DEMSim.DoDynamicsThenSync(0.);
                DEMSim.ChangeClumpFamily(1);
                size_t num_changed = 0;
                for (int i = 0; i < nW; i++) {
                    wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
                    float3 pos = ChVec2Float(wheel_pos[i]);
                    std::pair<float, float> Xrange =
                        std::pair<float, float>(pos.x - box_halfsize_x, pos.x + box_halfsize_x);
                    std::pair<float, float> Yrange =
                        std::pair<float, float>(pos.y - box_halfsize_y, pos.y + box_halfsize_y);
                    num_changed += DEMSim.ChangeClumpFamily(0, Xrange, Yrange);
                }
                std::cout << num_changed << " particles changed family number." << std::endl;
            }
        }

        // Run DEM
        d_start = std::chrono::high_resolution_clock::now();
        DEMSim.DoDynamics(step_size);
        d_end = std::chrono::high_resolution_clock::now();
        d_total += d_end - d_start;

        // Feed force
        if (curr_step % chrono_update_freq == 0) {
            for (int i = 0; i < nW; i++) {
                float3 F = trackers[i]->ContactAcc();
                F *= wheel_mass;
                float3 tor = trackers[i]->ContactAngAccLocal();
                tor = wheel_MOI * tor;

                wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
                Wheels[i]->Empty_forces_accumulators();
                Wheels[i]->Accumulate_force(Float2ChVec(F), wheel_pos[i], false);
                Wheels[i]->Accumulate_torque(Float2ChVec(tor), true);  // torque in DEME is local
                                                                       //
            }
            h_start = std::chrono::high_resolution_clock::now();
            sys.DoStepDynamics(step_size * chrono_update_freq);
            iris.Update();
            h_end = std::chrono::high_resolution_clock::now();
            h_total += h_end - h_start;
        }

        if (curr_step % report_steps == 0) {


            ChVector<float> rover_vel = Body_1->GetFrame_REF_to_abs().GetPos_dt();
            ChVector<float> rover_pos = Body_1->GetFrame_REF_to_abs().GetPos();
            ChQuaternion<float> rover_quat = Body_1->GetFrame_REF_to_abs().GetRot();


            max_v = max_v_finder->GetValue();
            std::cout << "Time: " << t << std::endl;
            std::cout << "rover pos: " << rover_pos.x() << ", " << rover_pos.y() << ", " << rover_pos.z() << std::endl;
            std::cout << "rover velo: " << rover_vel.x() << ", " << rover_vel.y() << ", " << rover_vel.z() << std::endl;
            std::cout << "rover quaternion: " << rover_quat.e0() << ", " << rover_quat.e1() << ", " << rover_quat.e2() << ", " << rover_quat.e3() << std::endl;
            std::cout << "motor torque: " << iris.GetDriveMotor(IrisWheelID::LF)->GetMotorTorque() << ", "  << iris.GetDriveMotor(IrisWheelID::RF)->GetMotorTorque() << ", " << iris.GetDriveMotor(IrisWheelID::LB)->GetMotorTorque() << ", " << iris.GetDriveMotor(IrisWheelID::RB)->GetMotorTorque()  << std::endl;

            for (int i = 0; i < nW; i++) {
                float3 F = trackers[i]->ContactAcc();
                F *= wheel_mass;
                float3 tor = trackers[i]->ContactAngAccLocal();
                tor = wheel_MOI * tor;


                std::cout << "wheel " << i << " f and trq from terrain: ";
                std::cout << F.x << ", " << F.y << ", " << F.z << ", " << tor.x << ", " << tor.y << ", " << tor.z << std::endl;            
            }

            std::cout << "========================" << std::endl;


            if (rover_pos.x() > world_x_size / 2. - 0.2 ) {
                std::cout << "This is far enough, stopping the simulation..." << std::endl;
                std::cout << "========================" << std::endl;
                DEMSim.DoDynamicsThenSync(0.);
                break;
            }
        }

    }
    std::cout << "Finishing up..." << std::endl;
    std::cout << h_total.count() << " seconds spent on host" << std::endl;
    std::cout << d_total.count() << " seconds spent on device" << std::endl;

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    return 0;
}




void PrepareParticles(DEMSolver& DEMSim, std::shared_ptr<DEMMaterial> mat_type_wheel){


    auto mat_type_wall = DEMSim.LoadMaterial({{"E", E}, {"nu", 0.3}, {"CoR", CoR}, {"mu", mu_wall}, {"Crr", 0.00}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", E}, {"nu", 0.3}, {"CoR", CoR}, {"mu", mu}, {"Crr", 0.00}});
    DEMSim.SetMaterialPropertyPair("mu", mat_type_wheel, mat_type_terrain, mu_wheel);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_wall, mat_type_terrain, mu_wall);

    DEMSim.InstructBoxDomainDimension(world_x_size, world_y_size, world_y_size);
    float bottom = -0.5;
    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_wall);
    DEMSim.AddBCPlane(make_float3(0,  world_y_size / 2, 0), make_float3(0, -1, 0), mat_type_wall);
    DEMSim.AddBCPlane(make_float3(0, -world_y_size / 2, 0), make_float3(0,  1, 0), mat_type_wall);
    // X-dir bounding planes
    DEMSim.AddBCPlane(make_float3(-world_x_size / 2., 0, 0), make_float3( 1, 0, 0), mat_type_wall);
    DEMSim.AddBCPlane(make_float3( world_x_size / 2., 0, 0), make_float3(-1, 0, 0), mat_type_wall);

    // Define the terrain particle templates
    // Calculate its mass and MOI
    float mass1 = 2.6e3 * 4.2520508;
    float3 MOI1 = make_float3(1.6850426, 1.6375114, 2.1187753) * 2.6e3;
    float mass2 = 2.6e3 * 2.1670011;
    float3 MOI2 = make_float3(0.57402126, 0.60616378, 0.92890173) * 2.6e3;
    // Scale the template we just created
    std::vector<double> scales = {0.0014, 0.00075833, 0.00044, 0.0003, 0.0002, 0.00018333, 0.00017};
    std::for_each(scales.begin(), scales.end(), [](double& r) { r *= 10.; });
    // Then load it to system
    std::shared_ptr<DEMClumpTemplate> my_template2 =
        DEMSim.LoadClumpType(mass2, MOI2, GetDEMEDataFile("clumps/triangular_flat_6comp.csv"), mat_type_terrain);
    std::shared_ptr<DEMClumpTemplate> my_template1 =
        DEMSim.LoadClumpType(mass1, MOI1, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {my_template2,
                                                                                DEMSim.Duplicate(my_template2),
                                                                                my_template1,
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1)};
    // Now scale those templates
    for (int i = 0; i < scales.size(); i++) {
        std::shared_ptr<DEMClumpTemplate>& my_template = ground_particle_templates.at(i);
        // Note the mass and MOI are also scaled in the process, automatically. But if you are not happy with this, you
        // can always manually change mass and MOI afterwards.
        my_template->Scale(scales.at(i));
        // Give these templates names, 0000, 0001 etc.
        char t_name[20];
        sprintf(t_name, "%04d", i);
        my_template->AssignName(std::string(t_name));
    }

    // Now we load clump locations from a checkpointed file
        std::cout << "Making terrain..." << std::endl;
        std::unordered_map<std::string, std::vector<float3>> clump_xyz;
        std::unordered_map<std::string, std::vector<float4>> clump_quaternion;
        try {
            clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_3e6.csv");
            clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_3e6.csv");
        } catch (...) {
            throw std::runtime_error("You will need to finish the GRCPrep demos first to obtain the checkpoint file GRC_3e6.csv, "
                             "in order to run this demo. That is a 4m by 2m GRC-1 terrain patch that is around 15cm "
                             "thick. If you don't have access to it, you can go to the forum "
                             "(https://groups.google.com/g/projectchrono) to ask the authors for it.");
        }
        std::vector<float3> in_xyz;   // particle position
        std::vector<float4> in_quat;  // particle quaternion 
        std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;  // particle clump type 
        unsigned int t_num = 0;
        for (int i = 0; i < scales.size(); i++) {
            char t_name[20];
            sprintf(t_name, "%04d", t_num);

            auto this_type_xyz = clump_xyz[std::string(t_name)];
            auto this_type_quat = clump_quaternion[std::string(t_name)];

            size_t n_clump_this_type = this_type_xyz.size();
            std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type
                      << std::endl;
            // Prepare clump type identification vector for loading into the system (don't forget type 0 in
            // ground_particle_templates is the template for rover wheel)
            std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                     ground_particle_templates.at(t_num));

            // Add them to the big long vector
            in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
            in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
            in_types.insert(in_types.end(), this_type.begin(), this_type.end());
            std::cout << "Added clump type " << t_num << std::endl;
            // Our template names are 0000, 0001 etc.
            t_num++;
        }

        // Finally, load the info into this batch
        DEMClumpBatch base_batch(in_xyz.size());
        base_batch.SetTypes(in_types);
        base_batch.SetPos(in_xyz);
        base_batch.SetOriQ(in_quat);

        DEMSim.AddClumps(base_batch);
}


//------------------------------------------------------------------
// Function to save the povray files of the MBD
//------------------------------------------------------------------
void SaveParaViewFiles(Iris& rover, path& rover_dir, unsigned int frame_number) {
    path filename;

    char f_name[20];
    sprintf(f_name, "%04d", frame_number);
    filename = rover_dir / ("./_" + std::string(f_name) + ".obj");

    std::vector<geometry::ChTriangleMeshConnected> meshes;

    // save the wheels to obj/vtk files
    for (int i = 0; i < 4; i++) {
        std::shared_ptr<ChBodyAuxRef> body;
        if (i == 0) {
            body = rover.GetWheel(IrisWheelID::LF)->GetBody();
        }
        if (i == 1) {
            body = rover.GetWheel(IrisWheelID::RF)->GetBody();
        }
        if (i == 2) {
            body = rover.GetWheel(IrisWheelID::LB)->GetBody();
        }
        if (i == 3) {
            body = rover.GetWheel(IrisWheelID::RB)->GetBody();
        }

        ChFrame<> body_ref_frame = body->GetFrame_REF_to_abs();
        ChVector<> body_pos = body_ref_frame.GetPos();      // body->GetPos();
        ChQuaternion<> body_rot = body_ref_frame.GetRot();  // body->GetRot();
        if (i == 1 || i == 3) {
            body_rot.Cross(body_rot, Q_from_AngZ(CH_C_PI));
        }

        auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        std::string obj_path = GetChronoDataFile("robot/iris/obj/iris_wheel.obj");
        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        meshes.push_back(*mmesh);
    }

    geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
}