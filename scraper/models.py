""" Database engine
"""

import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy import create_engine, exists, and_
from sqlalchemy import Column, Integer, Float, Text, Boolean, String, ForeignKey, UniqueConstraint

Base = declarative_base()
metadata = Base.metadata

db_file = 'sqlite:///ilthermo.db'

engine = create_engine(db_file, echo=False)
Session = sessionmaker(engine)
session = Session()


class Property(Base):
    __tablename__ = 'property'
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True)


class Paper(Base):
    __tablename__ = 'paper'
    id = Column(Integer, primary_key=True)
    year = Column(Integer)
    title = Column(Text, unique=True)
    author = Column(Text)
    journal = Column(Text)

    datas = relationship('Data', lazy='dynamic')


class UniqueIon(Base):
    __tablename__ = 'unique_ion'
    id = Column(Integer, primary_key=True)
    charge = Column(Integer)
    name = Column(Text, unique=True)
    smiles = Column(Text, unique=True)

    ions = relationship('Ion', lazy='dynamic')


class Ion(Base):
    __tablename__ = 'ion'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    searched = Column(Boolean)
    unique_ion_id = Column(Integer, ForeignKey(UniqueIon.id))

    unique_ion = relationship('UniqueIon', foreign_keys='Ion.unique_ion_id')


class UniqueMolecule(Base):
    __tablename__ = 'unique_molecule'
    __table_args__ = (UniqueConstraint('cation_id', 'anion_id', name='ion_id'),)
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    cation_id = Column(Integer, ForeignKey(UniqueIon.id))
    anion_id = Column(Integer, ForeignKey(UniqueIon.id))
    formula = Column(Text)
    smiles = Column(Text, unique=True)

    cation = relationship('UniqueIon', foreign_keys='UniqueMolecule.cation_id')
    anion = relationship('UniqueIon', foreign_keys='UniqueMolecule.anion_id')
    molecules = relationship('Molecule', lazy='dynamic')


class Molecule(Base):
    __tablename__ = 'molecule'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    code = Column(String(6))
    unique_molecule_id = Column(Integer, ForeignKey(UniqueMolecule.id))

    unique_molecule = relationship('UniqueMolecule', foreign_keys='Molecule.unique_molecule_id')


class Data(Base):
    __tablename__ = 'data'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(Molecule.id))
    paper_id = Column(Integer, ForeignKey(Paper.id))
    property_id = Column(Integer, ForeignKey(Property.id))
    phase = Column(String(20))
    t = Column(Float)
    p = Column(Float, nullable=True)
    value = Column(Float)
    stderr = Column(Float)

    molecule = relationship('Molecule', foreign_keys='Data.molecule_id')
    paper = relationship('Paper', foreign_keys='Data.paper_id')
    property = relationship('Property', foreign_keys='Data.property_id')


class DataSet(Base):
    __tablename__ = 'dataset'
    id = Column(Integer, primary_key=True)
    code = Column(String(5), unique=True)
    raw_data = Column(Text)


metadata.create_all(engine)
