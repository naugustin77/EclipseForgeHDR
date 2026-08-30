# Homebrew formula for a local install of EclipseForgeHDR.
#
# This installs from a source directory on disk rather than from a release
# tarball, so it is for people who have cloned the repo. Set ECLIPSEFORGE_SRC to
# the clone (or edit the url below), then:
#
#   ECLIPSEFORGE_SRC=/path/to/eclipseforgehdr \
#     brew install --formula ./formula/eclipseforgehdr.rb
#
# pipx is the simpler route and the one the README recommends; this exists for
# people who prefer to keep everything under brew.
#
class Eclipseforgehdr < Formula
  include Language::Python::Virtualenv

  desc "Corona HDR pipeline for total-solar-eclipse exposure brackets"
  homepage "https://github.com/USERNAME/eclipseforgehdr"
  url "file://#{ENV["ECLIPSEFORGE_SRC"] || File.expand_path("..", __dir__)}",
      using: :git, branch: "main"
  version "0.10.0"
  head "https://github.com/USERNAME/eclipseforgehdr.git", branch: "main"

  depends_on "python@3.12"

  def install
    virtualenv_install_with_resources
  end

  test do
    assert_match "eclipseforgehdr", shell_output("#{bin}/eclipseforgehdr --version")
  end
end
