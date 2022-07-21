### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 5004b766-9020-4f80-80e1-3a6469cc7a0c
begin
	using DataFrames
	using DelimitedFiles
	using StaticArrays
	using CSV
	using Interpolations
end

# ╔═╡ 4afea5e6-0860-11ed-1f7b-2999dc1d3750
md"""
# Imports
"""

# ╔═╡ 41abe6c8-f0db-4477-beb0-7e83973f4575
md"""
# File loading and parsing

Load all three coordinate files and parse data from parameter file. Sort final dataframe by patch number
"""

# ╔═╡ 2982cf60-155c-427f-b066-9c95f630f51f
begin
	sim_name = "run_cake_tests"
	root_directory = "/home/lucas/EinsteinToolkit/CactusX/exe/" * sim_name * "/"
	data_type = "vertex_coords"
	gf_name = "vcoord"
end;

# ╔═╡ 610e24ad-08e7-4bec-a615-6788ebb1afb5
begin
	files = @SVector[
		root_directory * "coordinates-" * data_type * ".it000000.x.tsv",
		root_directory * "coordinates-" * data_type * ".it000000.y.tsv",
		root_directory * "coordinates-" * data_type * ".it000000.z.tsv"
	]

	par_file = root_directory * sim_name * ".par"

	struct par_file_data
		outer_boundary_radius::Float64
		inner_boundary_radius::Float64
		cartesian_ncells::SVector{3, Int}
		angular_cells::Int
		radial_cells::Int
	end

	par_file_string = read(par_file, String)

	parameters = par_file_data(
		parse(Float64, match(r"MultiPatch::cake_outer_boundary_radius\s+=(\s+\d+.?\d+)", par_file_string)[1]),
		
		parse(Float64, match(r"MultiPatch::cake_inner_boundary_radius\s+=(\s+\d+.?\d+)", par_file_string)[1]),
		
		@SVector[
			parse(Int, match(r"MultiPatch::cake_cartesian_ncells_i\s+=\s+(\d+)", par_file_string)[1]),
			parse(Int, match(r"MultiPatch::cake_cartesian_ncells_j\s+=\s+(\d+)", par_file_string)[1]),
			parse(Int, match(r"MultiPatch::cake_cartesian_ncells_k\s+=\s+(\d+)", par_file_string)[1])
		],
		
		parse(Int, match(r"MultiPatch::cake_angular_cells\s+=\s+(\d+)", par_file_string)[1]),
		
		parse(Int, match(r"MultiPatch::cake_radial_cells\s+=\s+(\d+)", par_file_string)[1]),
	)

	coord_file_data = vcat(
		CSV.read(
			files[1],
			header=1,
			delim="\t",
			ignorerepeated=true,
			types = [Int, Float64, Int, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64],
			DataFrame
		),

		CSV.read(
			files[2],
			header=1,
			delim="\t",
			ignorerepeated=true,
			types = [Int, Float64, Int, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64],
			DataFrame
		),

		CSV.read(
			files[3],
			header=1,
			delim="\t",
			ignorerepeated=true,
			types = [Int, Float64, Int, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64],
			DataFrame
		),
	)
	coord_file_data = sort(coord_file_data, ["3:patch"])
end

# ╔═╡ 42a702d6-62b7-4fda-9e3f-68d332f32794
md"""
# Local to global coordinate transformations tests
Tests if `CarpetX` is doing the global to local coordinate transformations correctly. (No output is good)
"""

# ╔═╡ 7cfc0825-cc7a-40fb-a385-fd4288d6ee6b
function local2global(index::Int, a::Float64, b::Float64, c::Float64)
	r0 = parameters.inner_boundary_radius
	r1 = parameters.outer_boundary_radius

	if index == 0
		return @SVector [a, b, c]
	elseif index == 1
		return @SVector[
  			((1 - c)*r0 + (1 + c)*r1)/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
			(b*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
  			(a*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))
		]
	elseif index == 2
		@SVector[
  			-(((1 - c)*r0 + (1 + c)*r1)/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))),
  			-((b*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))),
  			(a*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))
		]
	elseif index == 3
		@SVector[
  			-((b*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))),
  			((1 - c)*r0 + (1 + c)*r1)/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
  			(a*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))
		]
	elseif index == 4
		@SVector[
  			(b*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
  			-(((1 - c)*r0 + (1 + c)*r1)/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))),
  			(a*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))
		]
	elseif index == 5
		@SVector[
  			-((a*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))),
  			(b*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
  			((1 - c)*r0 + (1 + c)*r1)/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c))
		]
	elseif index == 6
		@SVector[
  			(a*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
  			(b*((1 - c)*r0 + (1 + c)*r1))/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)),
  			-(((1 - c)*r0 + (1 + c)*r1)/sqrt(4 + 2*a^2*(1 + c) + 2*b^2*(1 + c)))
		]
	end
end

# ╔═╡ 417808e6-ef4b-4af7-8f43-3db3c749cc7c
for row in eachrow(coord_file_data)
	expected = local2global(
		row["3:patch"],
		row["8:x"],
		row["9:y"],
		row["10:z"],
	)
	
	computed = @SVector [
		row["11:" * gf_name * "x"],
		row["12:" * gf_name * "y"],
		row["13:" * gf_name * "z"]
	]

	if !(expected ≈ computed)
		println("Error")
	end
end

# ╔═╡ 91ffe50c-1eaf-431a-a93b-178ed5ab410e
md"""
# Global to local coordinate transformation tests
Tests if `CarpetX` is doinf the global to local coordinate transformations correctly (no output is good)
"""

# ╔═╡ f1e35e30-9239-4a81-945c-ec80b64aaab3
function global2local(index::Int, x::Float64, y::Float64, z::Float64)
	r0 = parameters.inner_boundary_radius
	r1 = parameters.outer_boundary_radius
	
	if index == 0
		return @SVector [x, y, z]
	elseif index == 1
		return @SVector[
  		z/x,
  		y/x,
  		(r0^2 - r1^2 + y^2 + z^2 + sqrt(4*r1^2*x^2 + (y^2 + z^2)^2 + 4*r0^2*(x^2 + y^2 + z^2) - 4*r0*r1*(2*x^2 + y^2 + z^2)))/(r0 - r1)^2
		]
	elseif index == 2
		return @SVector[
			-(z/x),
			y/x,
			(r0^2 - r1^2 + y^2 + z^2 + sqrt(4*r1^2*x^2 + (y^2 + z^2)^2 + 4*r0^2*(x^2 + y^2 + z^2) - 4*r0*r1*(2*x^2 + y^2 + z^2)))/(r0 - r1)^2
		]
	elseif index == 3
		return @SVector[
			z/y,
 		 	-(x/y),
  			(r0^2 - r1^2 + x^2 + z^2 + sqrt(4*r1^2*y^2 + (x^2 + z^2)^2 + 4*r0^2*(x^2 + y^2 + z^2) - 4*r0*r1*(x^2 + 2*y^2 + z^2)))/(r0 - r1)^2
		]
	elseif index == 4
		return @SVector[
			-(z/y),
  			-(x/y),
  			(r0^2 - r1^2 + x^2 + z^2 + sqrt(4*r1^2*y^2 + (x^2 + z^2)^2 + 4*r0^2*(x^2 + y^2 + z^2) - 4*r0*r1*(x^2 + 2*y^2 + z^2)))/(r0 - r1)^2
		]
	elseif index == 5
		return @SVector[
  			-(x/z),
  			y/z,
  			(r0^2 - r1^2 + x^2 + y^2 + sqrt((x^2 + y^2)*(4*r0*(r0 - r1) + x^2 + y^2) + 4*(r0 - r1)^2*z^2))/(r0 - r1)^2
		]
	elseif index == 6
		return @SVector[
			-(x/z),
  			-(y/z),
  			(r0^2 - r1^2 + x^2 + y^2 + sqrt((x^2 + y^2)*(4*r0*(r0 - r1) + x^2 + y^2) + 4*(r0 - r1)^2*z^2))/(r0 - r1)^2
		]
	end
end

# ╔═╡ 3abbdd40-8659-4488-ae53-ac58d5a7a8c4
for row in eachrow(coord_file_data)
	expected = global2local(
		row["3:patch"],
		row["11:" * gf_name * "x"],
		row["12:" * gf_name * "y"],
		row["13:" * gf_name * "z"]
	)
	
	computed = @SVector [
		row["8:x"],
		row["9:y"],
		row["10:z"],
	]

	if !(expected ≈ computed)
		println("Error")
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
Interpolations = "~0.14.0"
StaticArrays = "~1.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "00a19d6ab0cbdea2978fc23c5a6482e02c192501"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.0"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "e972716025466461a3dc1588d9168334b71aafff"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.1"

[[StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─4afea5e6-0860-11ed-1f7b-2999dc1d3750
# ╠═5004b766-9020-4f80-80e1-3a6469cc7a0c
# ╟─41abe6c8-f0db-4477-beb0-7e83973f4575
# ╠═2982cf60-155c-427f-b066-9c95f630f51f
# ╟─610e24ad-08e7-4bec-a615-6788ebb1afb5
# ╟─42a702d6-62b7-4fda-9e3f-68d332f32794
# ╟─7cfc0825-cc7a-40fb-a385-fd4288d6ee6b
# ╠═417808e6-ef4b-4af7-8f43-3db3c749cc7c
# ╟─91ffe50c-1eaf-431a-a93b-178ed5ab410e
# ╟─f1e35e30-9239-4a81-945c-ec80b64aaab3
# ╠═3abbdd40-8659-4488-ae53-ac58d5a7a8c4
# ╟─18744f01-8247-4b5e-819e-491176d931ae
# ╠═fea86325-ec57-4067-98d6-4e9b2e173f19
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
