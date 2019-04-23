# Microsoft Developer Studio Project File - Name="sheath" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=sheath - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "sheath.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "sheath.mak" CFG="sheath - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "sheath - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "sheath - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "sheath - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /stack:0x5f5e100 /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "sheath - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /stack:0xffffff /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "sheath - Win32 Release"
# Name "sheath - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\main.f90
DEP_F90_MAIN_=\
	".\Release\global.mod"\
	".\Release\grid_para.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod.f90
# End Source File
# Begin Source File

SOURCE=.\sub_boundary.f90
DEP_F90_SUB_B=\
	".\Release\constant.mod"\
	".\Release\diag.mod"\
	".\Release\global.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_dadi.f90
DEP_F90_SUB_D=\
	".\Release\constant.mod"\
	".\Release\discharge.mod"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\grid_para.mod"\
	".\Release\ion.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_diagnostic.f90
DEP_F90_SUB_DI=\
	".\Release\constant.mod"\
	".\Release\diag.mod"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\magnetic.mod"\
	".\Release\mcc_para.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_dielectric.f90
DEP_F90_SUB_DIE=\
	".\Release\constant.mod"\
	".\Release\global.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_DSMC.f90
DEP_F90_SUB_DS=\
	".\Release\constant.mod"\
	".\Release\DSMC.MOD"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\neutral.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_enter.f90
DEP_F90_SUB_E=\
	".\Release\discharge.mod"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_gather.f90
DEP_F90_SUB_G=\
	".\Release\constant.mod"\
	".\Release\diag.mod"\
	".\Release\global.mod"\
	".\Release\grid_para.mod"\
	".\Release\ion.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_gradient.f90
DEP_F90_SUB_GR=\
	".\Release\constant.mod"\
	".\Release\global.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\poisson.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_initial.f90
DEP_F90_SUB_I=\
	".\Release\constant.mod"\
	".\Release\diag.mod"\
	".\Release\discharge.mod"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\grid_para.mod"\
	".\Release\ion.mod"\
	".\Release\magnetic.mod"\
	".\Release\mcc_para.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\neutral.mod"\
	".\Release\poisson.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_move.f90
DEP_F90_SUB_M=\
	".\Release\constant.mod"\
	".\Release\diag.mod"\
	".\Release\global.mod"\
	".\Release\grid_para.mod"\
	".\Release\ion.mod"\
	".\Release\magnetic.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_output.f90
DEP_F90_SUB_O=\
	".\Release\constant.mod"\
	".\Release\diag.mod"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\grid_para.mod"\
	".\Release\ion.mod"\
	".\Release\magnetic.mod"\
	".\Release\mcc_para.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\neutral.mod"\
	".\Release\wall.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_smooth.f90
DEP_F90_SUB_S=\
	".\Release\grid.mod"\
	".\Release\mod_dadi.mod"\
	".\Release\mod_smooth.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sub_triangle_geometry.f90
# End Source File
# Begin Source File

SOURCE=.\sub_weighting.f90
DEP_F90_SUB_W=\
	".\Release\constant.mod"\
	".\Release\global.mod"\
	".\Release\grid.mod"\
	".\Release\grid_para.mod"\
	".\Release\magnetic.mod"\
	".\Release\neutral.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
