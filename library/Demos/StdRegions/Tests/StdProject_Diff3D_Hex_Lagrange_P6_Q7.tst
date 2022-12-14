<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff3D Hexahedron Lagrange basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s hexahedron -b GLL_Lagrange GLL_Lagrange GLL_Lagrange -o 6 6 6 -p 7 7 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.11314e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="4e-12">1.11413e-11</value>
        </metric>
    </metrics>
</test>


