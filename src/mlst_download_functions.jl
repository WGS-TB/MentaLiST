using Distributed
using LightXML
import GZip
@everywhere import GZip
using Printf

mentalist_shared_data_path = string(ENV["HOME"], "/", ".mentalist")

function _get_first_line(node)
  return split(content(node), "\n")[1]
end

function _older_than_a_day(file)
  return (time() - stat(file).ctime)/86400 > 1
end

function _pubmlst_xml(download_dir=mentalist_shared_data_path)
  if !isdir(download_dir)
    mkdir(download_dir)
  end
  dbases_file = string(download_dir, "/", "dbases.xml")
  if !isfile(dbases_file) || _older_than_a_day(dbases_file)
    @info("Downloading the MLST database xml file...")
    download("https://pubmlst.org/data/dbases.xml", dbases_file)
  end
  # get
  return LightXML.parse_file(dbases_file)
end

function _cgmlst_http(download_dir=mentalist_shared_data_path)
  if !isdir(download_dir)
    mkdir(download_dir)
  end
  cgmlst_file = string(download_dir, "/", "cgmlst.html")
  if !isfile(cgmlst_file)
    @info("Downloading the cgmlist HTML to find schema...")
    download("www.cgmlst.org/ncs", cgmlst_file)
  end
  return cgmlst_file
end


function list_pubmlst_schema(prefix)
  xroot = root(_pubmlst_xml())
  scheme_list = Tuple{Int,String}[]
  for (id,species) in enumerate(xroot["species"])
    sp_name = _get_first_line(species)
    if prefix == nothing || startswith(sp_name, prefix)
      push!(scheme_list,(id,sp_name))
    end
  end
  println("#id\torganism")
  for (id, sp_name) in scheme_list
    @printf "%d\t%-30s\n" id sp_name
  end
  @info("$(length(scheme_list)) scheme(s) found.")
end


@everywhere function _download_to_folder(url, output_dir, overwrite=false, filename=nothing)
  println("Trying to download: $url")
  filepath = joinpath(output_dir, filename == nothing ? basename(url) : filename)
  if overwrite || (!isfile(filepath) || _older_than_a_day(filepath))
    mkpath(output_dir)
    download(url, filepath)
  end
  return filepath
end

function _find_publmst_species(xml_xroot, target_id)
  for (id, species) in enumerate(xml_xroot["species"])
    sp_name = _get_first_line(species)
    if target_id == "$id" || target_id == sp_name
      return(species)
    end
  end
end

function download_pubmlst_scheme(target_species, output_dir, overwrite=false)
  loci_files = String[]
  xroot = root(_pubmlst_xml())
  @info("Searching for the scheme ... ")
  species = _find_publmst_species(xroot, target_species)
  if species == nothing
    @warn("I did not find this scheme on pubmlst, please check the species spelling or the ID and try again.")
    exit(-1)
    return
  end
  @info("Downloading scheme for $(_get_first_line(species)) ... ")
  db = species["mlst"][1]["database"][1]
  @info("Downloading profile ...")
  profile_url = content(db["profiles"][1]["url"][1])
  profile_file = _download_to_folder(profile_url, output_dir, overwrite)
  for locus in db["loci"][1]["locus"]
    locus_name = _get_first_line(locus)
    @info("Downloading locus $locus_name ...")
    locus_url = content(locus["url"][1])
    filepath = _download_to_folder(locus_url, output_dir)
    push!(loci_files, filepath)
  end
  @info("Finished downloading.")
  return (loci_files, profile_file)
end

function available_cgmlst_schema()
  scheme_list = Tuple{String, String}[]
  open(_cgmlst_http()) do f
    for l in eachline(f)
      for row in split(l, "<td>")
        m = match(r"<a href='https?://www.cgmlst.org/ncs/schema/(?<id>[0-9]+)/'><em>(?<species>.+)</em>", row)
        if m != nothing
          id, species = m.captures
          push!(scheme_list,(id,species))
        end
      end
    end
  end
  return scheme_list
end


function list_cgmlst_schema(prefix)
  println("#id\torganism")
  count = 0
  for (id, species) in available_cgmlst_schema()
    if prefix == nothing || startswith(species,prefix)
      count += 1
      @printf "%s\t%-30s\n" id species
    end
  end
  @info("$count schema found.")
end

function _find_cgmlst_id(target_id)
  for (id,species) in available_cgmlst_schema()
    if id == target_id || species == target_id
      return id
    end
  end
  return nothing
end

@everywhere function _download_enterobase_locus(locus, output_dir, url_folder)
  fasta_locus = joinpath(output_dir, "$locus.fa")
  if !isfile(fasta_locus)
    gzip_locus = _download_to_folder("https://enterobase.warwick.ac.uk/schemes/$url_folder/$locus.fasta.gz", output_dir, false, "$locus.fa.gz")
    # gunzip to a FASTA and remove the gzip file;
    f_in = GZip.open(gzip_locus)
    f_out = open(fasta_locus, "w")
    while !eof(f_in)
      write(f_out, readline(f_in))
    end
    close(f_in)
    close(f_out)
    rm(gzip_locus)
  end
  return fasta_locus
end

function download_enterobase_scheme(scheme, s_type, output_dir, overwrite=false)
  verbose = Dict("S"=>"Salmonella", "Y"=>"Yersinia", "E"=>"Escherichia/Shigella", "C" => "clostridium")
  sp = Dict("E"=>"Escherichia", "S"=>"Salmonella", "Y"=>"Yersinia", "C" => "clostridium")

  if !haskey(sp, scheme)
    @info("Scheme has to be E, S, Y or C.")
    exit(-1)
  end
  # Version: empty for wg, v2 for Salmonella cg, v1 for other cg.
  version = s_type == "wg" ? "" : (scheme == "S" ? "v2" : "v1")
  scheme_name = "$(s_type)MLST$(version)"
  url_folder = "$(sp[scheme]).$(scheme_name)"
  filename = joinpath(dirname(@__FILE__), "../scripts/$(sp[scheme]).txt")

  loci = String[]
  open(filename) do f
    for l in eachline(f)
      values = split(strip(l),"\t")
      if startswith(values[1], s_type) # add locus if it is cg or wgmlst, according to input option
        push!(loci, split(values[2],' ')[1])
      end
    end
  end
  loci_files = pmap(locus->_download_enterobase_locus(locus, output_dir, url_folder), loci)
  return loci_files
end

function download_cgmlst_scheme(target_id, output_dir)
  id = _find_cgmlst_id(target_id)
  if id == nothing
    @warn("Id/species ($target_id) not found!")
    exit(-1)
  end
  @info("Downloading cgMLST scheme ...")
  scheme_zip_file = _download_to_folder("https://www.cgmlst.org/ncs/schema/$id/alleles", output_dir)
  locus_files = String[]
  locus = ""
  @info("Unzipping cgMLST scheme into individual FASTA files for each locus ...")
  scheme_dirname = dirname(scheme_zip_file)
  run(`unzip -oq $scheme_zip_file -d $scheme_dirname/tmp`)
  rm(scheme_zip_file)
  scheme_files = readdir(joinpath(scheme_dirname, "tmp"))

  for scheme_file in scheme_files
    # print to show progress:
    if length(locus_files) % 200 == 0
        print(".")
    end
    scheme_file_path = joinpath(scheme_dirname, scheme_file)
    push!(locus_files, scheme_file_path)
    # get locus ID from filename
    locus = split(scheme_file, ".")[1]
    fh = open(scheme_file_path, "w")
    for l in eachline(joinpath(scheme_dirname, "tmp", scheme_file))
      if length(strip(l)) == 0
        continue
      end
      if l[1] == '>'
        l = '>' * locus * "_" * l[2:end] * "\n"
        write(fh, l)
      else
        write(fh, l * "\n")
      end
    end
    close(fh)
  end
  println()
  total_loci = length(locus_files)
  @info("$total_loci loci found.")
  rm(joinpath(scheme_dirname, "tmp"), recursive=true)
  return locus_files
end
