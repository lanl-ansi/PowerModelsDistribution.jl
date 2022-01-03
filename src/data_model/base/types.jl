
abstract type DataModel end

# Core Data Model Types
abstract type OpenDssRawModel <: DataModel end
abstract type OpenDssModel <: DataModel end
abstract type EngineeringModel <: DataModel end
abstract type MathematicalModel <: DataModel end

# OpenDss Model Abstract Object Types
abstract type OpenDssObject <: OpenDssModel end
abstract type OpenDssNodeObject <: OpenDssObject end
abstract type OpenDssEdgeObject <: OpenDssObject end
abstract type OpenDssDataObject <: OpenDssObject end
abstract type OpenDssControlObject <: OpenDssObject end

# Engineering Model Abstract Object Types
abstract type EngObject <: EngineeringModel end
abstract type EngNodeObject <: EngObject end
abstract type EngEdgeObject <: EngObject end
abstract type EngDataObject <: EngObject end
abstract type EngControlObject <: EngObject end

# Mathematical Model Abstract Object Types
abstract type MathObject <: MathematicalModel end
abstract type MathNodeObject <: MathObject end
abstract type MathEdgeObject <: MathObject end
abstract type MathDataObject <: MathObject end
abstract type MathControlObject <: MathObject end
