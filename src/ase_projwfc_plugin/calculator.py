"""ASE calculator for PW2Wannier90."""

import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

from ase import Atoms
from ase.calculators.genericfileio import BaseProfile, CalculatorTemplate, GenericFileIOCalculator
from ase.io import read, write


class PW2Wannier90Profile(BaseProfile):
    """Profile for PW2Wannier90 calculator."""

    def get_calculator_command(self, inputfile: str) -> list[str]:
        """Return the arguments to run the calculator."""
        return [inputfile]

    def version(self) -> str:
        """Return the version of the PW2Wannier90 executable."""
        raise NotImplementedError()


class PW2Wannier90Template(CalculatorTemplate):
    """Template for PW2Wannier90 calculator."""

    _label = "projwfc"

    def __init__(self) -> None:
        super().__init__("projwfc", ["spreads", "band structure"])
        self.inputname = self._label + ".p2wi"
        self.outputname = self._label + ".p2wo"
        self.errorname = self._label + ".p2we"

    def write_input(
        self,
        profile: PW2Wannier90Profile,
        directory: Path,
        atoms: Atoms,
        parameters: dict[str, Any],
        properties: list[str],
    ) -> None:
        """Generate a PW2Wannier90 input file."""
        dst = directory / self.inputname

        write(dst, atoms, format="projwfc-in", **parameters)

    def execute(self, directory: Path, profile: PW2Wannier90Profile) -> None:
        """Run the PW2Wannier90 calculation."""
        profile.run(directory, self.inputname, self.outputname, errorfile=self.errorname)

    def read_results(self, directory: os.PathLike[Any]) -> Mapping[str, Any]:
        """Read results from a PW2Wannier90 calculation."""
        directory_path = Path(directory)
        atoms_in = read(directory_path / self.inputname, format="projwfc-in")
        if not isinstance(atoms_in, Atoms):
            raise ValueError("Input file did not return an Atoms object")
        
        atoms_out = read(directory_path / self.outputname, format="projwfc-out")
        if not isinstance(atoms_out, Atoms):
            raise ValueError("Output file did not return an Atoms object")
        if atoms_out.calc is None:
            raise ValueError("Output Atoms object has no calculator attached")

        results = atoms_out.calc.results
        if not isinstance(results, dict):
            raise ValueError("Output file did not return a results dictionary")
    
        if atoms_in.calc.parameters.get("bands_plot", False):
            bs = _read_bandstructure(
                directory_path, self._label, atoms_in.cell, atoms_in.calc.parameters["kpoint_path"]
            )
            if bs is not None:
                results["band_structure"] = bs

        return results

    def load_profile(self, cfg: dict[str, Any], **kwargs: Any) -> PW2Wannier90Profile:
        """Load a PW2Wannier90Profile from a configuration."""
        return cast(PW2Wannier90Profile, PW2Wannier90Profile.from_config(cfg, self.name, **kwargs))


class PW2Wannier90(GenericFileIOCalculator):
    """An ASE calculator for PW2Wannier90."""

    def __init__(
        self, profile: PW2Wannier90Profile | None = None, directory: str = ".", **kwargs: Any
    ):
        template = PW2Wannier90Template()
        super().__init__(profile=profile, template=template, directory=directory, parameters=kwargs)

