##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=Deduplifer
ConfigurationName      :=Release
WorkspacePath          :=/home/arunachala/Documents/BioSrc
ProjectPath            :=/home/arunachala/Documents/BioSrc/Deduplifer
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Oliver Rimington
Date                   :=17/08/16
CodeLitePath           :=/home/arunachala/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=$(PreprocessorSwitch)NDEBUG 
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="Deduplifer.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -O2 -Wall $(Preprocessors)
CFLAGS   :=  -O2 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/Utility.cpp$(ObjectSuffix) $(IntermediateDirectory)/Scaffold.cpp$(ObjectSuffix) $(IntermediateDirectory)/DigestBlast.cpp$(ObjectSuffix) $(IntermediateDirectory)/ScaffNet.cpp$(ObjectSuffix) $(IntermediateDirectory)/st_ed.cpp$(ObjectSuffix) $(IntermediateDirectory)/edge_link.cpp$(ObjectSuffix) $(IntermediateDirectory)/DigestGenome.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Release || $(MakeDirCommand) ./Release


$(IntermediateDirectory)/.d:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM "main.cpp"

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) "main.cpp"

$(IntermediateDirectory)/Utility.cpp$(ObjectSuffix): Utility.cpp $(IntermediateDirectory)/Utility.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/Utility.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Utility.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Utility.cpp$(DependSuffix): Utility.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Utility.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Utility.cpp$(DependSuffix) -MM "Utility.cpp"

$(IntermediateDirectory)/Utility.cpp$(PreprocessSuffix): Utility.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Utility.cpp$(PreprocessSuffix) "Utility.cpp"

$(IntermediateDirectory)/Scaffold.cpp$(ObjectSuffix): Scaffold.cpp $(IntermediateDirectory)/Scaffold.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/Scaffold.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Scaffold.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Scaffold.cpp$(DependSuffix): Scaffold.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Scaffold.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Scaffold.cpp$(DependSuffix) -MM "Scaffold.cpp"

$(IntermediateDirectory)/Scaffold.cpp$(PreprocessSuffix): Scaffold.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Scaffold.cpp$(PreprocessSuffix) "Scaffold.cpp"

$(IntermediateDirectory)/DigestBlast.cpp$(ObjectSuffix): DigestBlast.cpp $(IntermediateDirectory)/DigestBlast.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/DigestBlast.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/DigestBlast.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/DigestBlast.cpp$(DependSuffix): DigestBlast.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/DigestBlast.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/DigestBlast.cpp$(DependSuffix) -MM "DigestBlast.cpp"

$(IntermediateDirectory)/DigestBlast.cpp$(PreprocessSuffix): DigestBlast.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/DigestBlast.cpp$(PreprocessSuffix) "DigestBlast.cpp"

$(IntermediateDirectory)/ScaffNet.cpp$(ObjectSuffix): ScaffNet.cpp $(IntermediateDirectory)/ScaffNet.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/ScaffNet.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/ScaffNet.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/ScaffNet.cpp$(DependSuffix): ScaffNet.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/ScaffNet.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/ScaffNet.cpp$(DependSuffix) -MM "ScaffNet.cpp"

$(IntermediateDirectory)/ScaffNet.cpp$(PreprocessSuffix): ScaffNet.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/ScaffNet.cpp$(PreprocessSuffix) "ScaffNet.cpp"

$(IntermediateDirectory)/st_ed.cpp$(ObjectSuffix): st_ed.cpp $(IntermediateDirectory)/st_ed.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/st_ed.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/st_ed.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/st_ed.cpp$(DependSuffix): st_ed.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/st_ed.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/st_ed.cpp$(DependSuffix) -MM "st_ed.cpp"

$(IntermediateDirectory)/st_ed.cpp$(PreprocessSuffix): st_ed.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/st_ed.cpp$(PreprocessSuffix) "st_ed.cpp"

$(IntermediateDirectory)/edge_link.cpp$(ObjectSuffix): edge_link.cpp $(IntermediateDirectory)/edge_link.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/edge_link.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/edge_link.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/edge_link.cpp$(DependSuffix): edge_link.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/edge_link.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/edge_link.cpp$(DependSuffix) -MM "edge_link.cpp"

$(IntermediateDirectory)/edge_link.cpp$(PreprocessSuffix): edge_link.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/edge_link.cpp$(PreprocessSuffix) "edge_link.cpp"

$(IntermediateDirectory)/DigestGenome.cpp$(ObjectSuffix): DigestGenome.cpp $(IntermediateDirectory)/DigestGenome.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/arunachala/Documents/BioSrc/Deduplifer/DigestGenome.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/DigestGenome.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/DigestGenome.cpp$(DependSuffix): DigestGenome.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/DigestGenome.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/DigestGenome.cpp$(DependSuffix) -MM "DigestGenome.cpp"

$(IntermediateDirectory)/DigestGenome.cpp$(PreprocessSuffix): DigestGenome.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/DigestGenome.cpp$(PreprocessSuffix) "DigestGenome.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Release/


