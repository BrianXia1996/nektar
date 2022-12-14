<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_IP_GAUSS.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_IP_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-5">5.50779e+00</value>
            <value variable="rhou" tolerance="1e-12">2.01632e+03</value>
            <value variable="rhov" tolerance="5e-4">2.47024e+01</value>
            <value variable="E" tolerance="1e-12">6.27023e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-5">3.2956e-01</value>
            <value variable="rhou" tolerance="2e-4">8.89524e+01</value>
            <value variable="rhov" tolerance="1e-8">1.84552e+01</value>
            <value variable="E" tolerance="1e-12">2.87040e+05</value>
        </metric>
    </metrics>
</test>


