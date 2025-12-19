from pydantic import BaseModel, Field, EmailStr, ConfigDict
from typing import Optional, List, Dict, Any
from bson import ObjectId
from datetime import datetime

# -----------------------------------------------------
# Helper: Pydantic v2 compatible ObjectId wrapper
# -----------------------------------------------------
class PyObjectId(ObjectId):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        if isinstance(v, ObjectId):
            return v
        if isinstance(v, str):
            try:
                return ObjectId(v)
            except Exception:
                raise ValueError("Invalid ObjectId")
        raise TypeError("ObjectId must be a string or ObjectId")

    @classmethod
    def __get_pydantic_json_schema__(cls, core_schema, handler):
        json_schema = handler(core_schema)
        json_schema.update(type="string")
        return json_schema


# -----------------------------------------------------
# Base Mongo Document
# -----------------------------------------------------
class MongoBaseModel(BaseModel):
    # ❗ No "_id" defined here
    model_config = ConfigDict(
        populate_by_name=True,
        arbitrary_types_allowed=True,
        json_encoders={ObjectId: str},
    )



# -----------------------------------------------------
# USER MODELS
# -----------------------------------------------------
class UserCreate(BaseModel):
    username: str
    email: EmailStr
    password: str


class UserDB(MongoBaseModel):
    username: str
    email: EmailStr
    hashed_password: str
    created_at: datetime = Field(default_factory=datetime.utcnow)


class UserResponse(BaseModel):
    id: str
    username: str
    email: EmailStr


# -----------------------------------------------------
# PROJECT MODELS
# -----------------------------------------------------
class ProjectCreate(BaseModel):
    name: str
    description: Optional[str] = None


class ProjectDB(MongoBaseModel):
    id: PyObjectId = Field(default_factory=PyObjectId, alias="_id")
    user_id: str
    name: str
    description: Optional[str] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)



class ProjectResponse(BaseModel):
    id: str
    name: str
    description: Optional[str] = None
    created_at: datetime


# -----------------------------------------------------
# MOLECULE MODELS
# -----------------------------------------------------
class MolProperties(BaseModel):
    logp: float
    qed: float
    mw: float
    tpsa: float
    hbd: int
    hba: int
    rot_bonds: int


class AdmetProperties(BaseModel):
    absorption: float
    distribution: float
    metabolism: float
    excretion: float
    toxicity: float


class MoleculeCreate(BaseModel):
    project_id: str
    name: str
    smiles: str
    generated_by: str
    tags: List[str] = []
    properties: MolProperties
    admet: Optional[AdmetProperties] = None
    tox_alerts: Optional[Dict[str, Any]] = {}


class MoleculeDB(MongoBaseModel):
    user_id: str
    project_id: str
    name: str
    smiles: str
    generated_by: str
    tags: List[str]
    properties: MolProperties
    admet: Optional[AdmetProperties] = None
    tox_alerts: Optional[Dict[str, Any]] = {}
    created_at: datetime = Field(default_factory=datetime.utcnow)


class MoleculeResponse(BaseModel):
    id: str
    project_id: str
    name: str
    smiles: str
    generated_by: str
    tags: List[str]
    properties: MolProperties
    admet: Optional[AdmetProperties]
    tox_alerts: Optional[Dict[str, Any]]
    created_at: datetime


# -----------------------------------------------------
# TOKEN MODELS
# -----------------------------------------------------
class Token(BaseModel):
    access_token: str
    refresh_token: str
    token_type: str = "bearer"


class TokenPayload(BaseModel):
    sub: str
    exp: int

# -----------------------------------------------------
# JOB MODELS
# -----------------------------------------------------
class JobCreate(BaseModel):
    task_type: str  # one of: GENERATE_RANDOM, GENERATE_TARGETED, OPTIMIZE_LEAD
    params: Dict[str, Any]

class JobDB(MongoBaseModel):
    id: PyObjectId = Field(default_factory=PyObjectId, alias="_id")
    user_id: str
    task_type: str
    status: str = "PENDING"  # PENDING, PROCESSING, COMPLETED, FAILED
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    params: Dict[str, Any]
    result: Optional[List[Any]] = None
    error: Optional[str] = None

class JobResponse(BaseModel):
    job_id: str
    status: str
    created_at: datetime
    updated_at: datetime
    result: Optional[List[Any]] = None
    error: Optional[str] = None
