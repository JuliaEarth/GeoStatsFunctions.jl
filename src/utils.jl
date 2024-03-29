# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    spheredir(θ, φ)

Return the 3D direction given polar angle `θ` and
azimuthal angle `φ` in degrees according to the ISO
convention.
"""
function spheredir(theta, phi)
  θ, φ = deg2rad(theta), deg2rad(phi)
  Vec(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
end
