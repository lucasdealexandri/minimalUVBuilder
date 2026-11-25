(* ::Package:: *)
(* Dump the Mathematica and FeynRules symbol lists used by tools/fr_names.py as reserved names.

   Writes  src/uvbuilder/tools/reserved/mathematica_<version>.txt  and  src/uvbuilder/tools/reserved/feynrules_<version>.txt
   (one symbol per line, '#' header lines with versions and date). All files present are unioned
   by tools/fr_names.language_lists, so lists from several versions can coexist.

   Command line (from the repository root):
       wolframscript -file scripts/dump_reserved_symbols.wl /path/to/feynrules-current
   Notebook: set frPath and outDir below by hand and evaluate the file's cells. *)

frPath = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], Environment["FEYNRULES_PATH"]];
outDir = If[StringQ[$InputFileName] && $InputFileName =!= "",
   FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "uvbuilder", "tools", "reserved"}],
   FileNameJoin[{Directory[], "src", "uvbuilder", "tools", "reserved"}]];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

header[what_] := {"# " <> what, "# Mathematica " <> $Version, "# generated " <> DateString[{"Year", "-", "Month", "-", "Day"}]};

Export[FileNameJoin[{outDir, "mathematica_" <> ToString[$VersionNumber] <> ".txt"}],
  Join[header["Names[\"System`*\"]"], Names["System`*"]], "Lines"];

If[!StringQ[frPath] || !DirectoryQ[frPath],
  Print["FeynRules path not given (argument or FEYNRULES_PATH): only the Mathematica list was written."],
  $FeynRulesPath = SetDirectory[frPath];
  << FeynRules`;
  (* Resolved by name after loading: a literal symbol here would be created in Global` when this
     If[...] is parsed, i.e. before FeynRules defines it. *)
  frVersion = ToString[ToExpression["FR$VersionNumber"]];
  If[!StringMatchQ[frVersion, RegularExpression["[0-9][0-9.]*"]], frVersion = "unknown"];
  Export[FileNameJoin[{outDir, "feynrules_" <> frVersion <> ".txt"}],
    Join[header["Names[\"FeynRules`*\"]"], {"# FeynRules " <> frVersion}, Names["FeynRules`*"]], "Lines"]
];
Print["Reserved-symbol lists written to ", outDir];
