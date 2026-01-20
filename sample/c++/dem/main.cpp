#include <particle_simulator.hpp>
#include <random>

#include "force.h"
#include "class.h"
#include "surface.h"
//#include "init/InfiniteCylinder2.h"
//#include "init/InfiniteCylinder2_vibrating.h"
#include "init/RectangularBox_vibrating.h"
#include "integral.h"
#include "io.h"

void RemoveParticles(PS::ParticleSystem<FP>& ptcl) {
    // 距離しきい値: 壁からこの距離以上離れていれば削除
    constexpr PS::F64 kRemoveDistance = 0.002;  // 2 mm
    std::vector<PS::S32> remove_ids;
    remove_ids.reserve(ptcl.getNumberOfParticleLocal());
    for (int i = 0; i < ptcl.getNumberOfParticleLocal(); ++i) {
        const PS::F64 ax = std::abs(ptcl[i].pos.x) - kBoxHalfX;
        const PS::F64 ay = std::abs(ptcl[i].pos.y) - kBoxHalfY;
        const PS::F64 az = std::abs(ptcl[i].pos.z) - kBoxHalfZ;
        const PS::F64 dx = std::max(ax, 0.0);
        const PS::F64 dy = std::max(ay, 0.0);
        const PS::F64 dz = std::max(az, 0.0);
        // 直方体外側の最短距離（外側にいるときのみ正）
        const PS::F64 dist_out = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (dist_out >= kRemoveDistance) {
            remove_ids.push_back(static_cast<PS::S32>(ptcl[i].id));
        }
    }

    if (!remove_ids.empty()) {
        ptcl.removeParticle(&remove_ids[0], remove_ids.size());
    }
}

int main(int argc, char* argv[]) {
    PS::Initialize(argc, argv);
    PS::ParticleSystem<FP> ptcl;
    ptcl.initialize();
    PS::TreeForForceShort<Force, FP, FP>::Symmetry tree_coll;

    PS::DomainInfo dinfo;
    dinfo.initialize();
    system_t sysinfo;
    FileIO<FP> io;

    if (argc == 1) {
        Problem<FP>::SetIC(ptcl, dinfo, sysinfo);
        tree_coll.initialize(ptcl.getNumberOfParticleLocal(), 0.5, 1, 1);
    } else {
        io.Restore(ptcl, sysinfo);
        tree_coll.initialize(ptcl.getNumberOfParticleLocal(), 0.5, 1, 1);
        sysinfo.end_time = 10.0;
        goto savepoint;
    }

    // for (int i = 0; i < ptcl.getNumberOfParticleLocal(); i++) {
    //     std::cout << "i= " << i << " pos= " << ptcl[i].pos << " rad= " << ptcl[i].rad << " kind= " << ptcl[i].kind << std::endl;
    // }

    dinfo.decomposeDomainAll(ptcl);
    ptcl.exchangeParticle(dinfo);
    tree_coll.calcForceAllAndWriteBack(Force(), ptcl, dinfo);
    Problem<FP>::externalForce(ptcl, sysinfo);

    for (sysinfo.time = 0, sysinfo.step = 0; sysinfo.time < sysinfo.end_time; sysinfo.time += sysinfo.dt, ++sysinfo.step) {
        dinfo.decomposeDomainAll(ptcl);
        ptcl.exchangeParticle(dinfo);
        sysinfo.dt = GetGlobalTimestep<FP>(ptcl);
        for (int i = 0; i < ptcl.getNumberOfParticleLocal(); ++i) {
            ptcl[i].kick(sysinfo.dt);
            ptcl[i].drift(sysinfo.dt);
            ptcl[i].clear();
        }
        tree_coll.calcForceAllAndWriteBack(Force(), ptcl, dinfo);
        Problem<FP>::externalForce(ptcl, sysinfo);
        for (int i = 0; i < ptcl.getNumberOfParticleLocal(); ++i) {
            ptcl[i].kick2(sysinfo.dt);
        }
        Problem<FP>::postTimestep(ptcl, sysinfo);
        ptcl.adjustPositionIntoRootDomain(dinfo);
        if (PS::Comm::getRank() == 0 && sysinfo.step % 100 == 0)
            std::cerr << "time = " << sysinfo.time << " (dt = " << sysinfo.dt << ")" << std::endl;
        io.OutputFileWithTimeInterval(ptcl, sysinfo);
        if (sysinfo.step % 10000 == 0) io.Create(ptcl, sysinfo);
    savepoint:;
        RemoveParticles(ptcl);
    }

    PS::Finalize();
    return 0;
}
