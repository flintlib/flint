{
  sources ? import ./npins,
  pkgs ? import sources.nixpkgs { },
}:

let
  inherit (pkgs.lib) genAttrs;

  packages = with pkgs; [
    # Dependencies
    gmp
    mpfr
    ntl
    openblas

    # Build tools
    autoconf
    automake
    gnumake
    libtool
  ];

  stdenvs = [
    "clang"
    "gcc9"
    "gcc10"
    "gcc11"
    "gcc12"
    "gcc13"
    "gcc14"
  ];
in

{
  devShell = pkgs.mkShell {
    name = "flint.dev";
    inherit packages;
    passthru = genAttrs stdenvs (
      stdenv:
      (pkgs.mkShell.override { stdenv = pkgs.${stdenv + "Stdenv"}; }) {
        name = "flint.dev-${stdenv}";
        inherit packages;
      }
    );
  };
}
