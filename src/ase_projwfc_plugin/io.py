"""Reads projwfc files.

"""

import json
from typing import Any, ClassVar, Generator, TextIO
from ase import Atoms
from ase.calculators.calculator import Calculator, PropertyNotImplementedError
from ase.utils.plugins import ExternalIOFormat
from pydantic_espresso.models.qe_7_5.projwfc import PW2WANNIER90EspressoInput as PW2Wannier90Input
from pint import Quantity

from ase_projwfc_plugin.utils import unit_registry

class PW2Wannier90ResultsCalculator(Calculator):
    """Container for PW2Wannier90 results."""

    implemented_properties: ClassVar[list[str]] = []  # type: ignore[misc]
    # Ignoring type override because ASE wrongly does not declare this as a ClassVar

    def __init__(self, atoms: Atoms, **results: Any):
        super().__init__(atoms=atoms.copy())
        self.results: dict[str, Any] = {}
        for key, value in results.items():
            if value in [None, []]:
                continue
            self.results[key] = value

    def get_property(
        self, name: str, atoms: Atoms | None = None, allow_calculation: bool = True
    ) -> Any:
        """Return a property."""
        if atoms is None:
            atoms = self.atoms
        if name not in self.results or self.check_state(atoms):
            if allow_calculation:
                raise PropertyNotImplementedError(f"The property {name} is not available.")
            return None
        return self.results[name] 

def _line_to_keyval(line: str) -> tuple[str, Any]:
    """Convert a line in an input file to a key-value pair."""
    if '=' not in line:
        raise ValueError(f'Line "{line}" is not a valid key-value pair')
    key, val = line.split('=', 1)
    key = key.strip().lower()
    val = val.strip('"').strip("'")
    try:
        # Try to convert to int/float/bool as appropriate
        return key, json.loads(val)
    except json.JSONDecodeError:
        return key, val

def read_projwfc_in(fileobj: str | TextIO) -> Atoms:
    """Parse a projwfc input file

    inputs are a fortran-namelist format with custom blocks of data.
    The namelist is parsed as a dict and an atoms object is constructed
    from the included information.

    Parameters
    ----------
    fileobj : file | str
        A file-like object that supports line iteration with the contents
        of the input file, or a filename.

    Returns
    -------
    atoms : Atoms
        Structure defined in the input file.

    Raises
    ------
    KeyError
        Raised for missing keys that are required to process the file
    """
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'rU')

    flines = fileobj.readlines()
    if flines[0].strip().lower() != '&inputpp':
        raise ValueError('Not a valid projwfc input file')
    if flines[-1].strip() != '/':
        raise ValueError('projwfc input file not terminated correctly')
    
    data = {}
    for line in flines[1:-1]:
        key, value = _line_to_keyval(line)
        data[key] = value

    atoms = Atoms()
    calc = PW2Wannier90ResultsCalculator(atoms=atoms)
    calc.parameters.update(**data)
    atoms.calc = calc

    return atoms


def write_projwfc_in(fd: str | TextIO, atoms: Atoms, **kwargs: Any) -> None:
    """
    Create an input file for projwfc.

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.

    """

    model = PW2Wannier90Input(inputpp=kwargs)
    fd.write(str(model))


def read_projwfc_out(fd: str | TextIO) -> Generator[Atoms, None, None]:
    """
    Reads projwfc output files

    Parameters
    ----------
    fd : file|str
        A file like object or filename

    Yields
    ------
    structure : atoms
        An Atoms object with an attached SinglePointCalculator containing
        any parsed results
    """

    if isinstance(fd, str):
        fd = open(fd, 'rU')

    flines = fd.readlines()

    structure = Atoms()

    job_done = False
    walltime: Quantity | None = None
    amn_calculated = False
    mmn_calculated = False

    for line in flines:
        if 'JOB DONE' in line:
            job_done = True
        elif line.strip().startswith('PW2WANNIER'):
            time_str = line.split()[-2]
            walltime = unit_registry(time_str)
        elif 'AMN calculated' in line:
            amn_calculated = True
        elif 'MMN calculated' in line:
            mmn_calculated = True

    calc = PW2Wannier90ResultsCalculator(atoms=structure)
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    calc.results['amn_calculated'] = amn_calculated
    calc.results['mmn_calculated'] = mmn_calculated

    structure.calc = calc

    yield structure

PW2WANNIER90_INPUT_FORMAT = ExternalIOFormat(
    "Input file for projwfc", "1F", module="ase_projwfc_plugin.io", ext=".p2wi"
)

PW2WANNIER90_OUTPUT_FORMAT = ExternalIOFormat(
    "Output file for projwfc", "1F", module="ase_projwfc_plugin.io", ext=".p2wo"
)