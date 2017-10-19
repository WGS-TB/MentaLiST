using Suppressor
@suppress_err begin
  using LightXML
  import GZip
end

function _get_first_line(node)
  return split(content(node), "\n")[1]
end

function _older_than_a_day(file)
  return (time() - stat(file).ctime)/86400 > 1
end
function _pubmlst_xml()
  if !isfile("dbases.xml") || _older_than_a_day("dbases.xml")
    info("Downloading the MLST database xml file...")
    download("https://pubmlst.org/data/dbases.xml", "dbases.xml")
  end
  # get
  return LightXML.parse_file("dbases.xml")
end

function _cgmlst_http()
  if !isfile("cgmlst.html")
    info("Downloading the cgmlist HTML to find schema...")
    download("www.cgmlst.org/ncs", "cgmlst.html")
  end
  return "cgmlst.html"
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
  for (id, sp_name) in scheme_list
    @printf "%-30s ID:%d\n" sp_name id
  end
  println("$(length(scheme_list)) schema found.")
end


function _download_to_folder(url, output_dir, overwrite=false)
  filepath = joinpath(output_dir, basename(url))
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
  info("Searching for the scheme ... ")
  species = _find_publmst_species(xroot, target_species)
  if species == nothing
    Lumberjack.warn("I did not found this scheme on pubmlst, please check the species spelling or the ID and try again.")
    exit(-1)
    return
  end
  info("Downloading scheme for $(_get_first_line(species)) ... ")
  db = species["mlst"][1]["database"][1]
  info("Downloading profile ...")
  profile_url = content(db["profiles"][1]["url"][1])
  profile_file = _download_to_folder(profile_url, output_dir, overwrite)
  for locus in db["loci"][1]["locus"]
    locus_name = _get_first_line(locus)
    info("Downloading locus $locus_name ...")
    locus_url = content(locus["url"][1])
    filepath = _download_to_folder(locus_url, output_dir)
    push!(loci_files, filepath)
  end
  info("Finished downloading.")
  return (loci_files, profile_file)
end

function available_cgmlst_schema()
  scheme_list = Tuple{String, String}[]
  open(_cgmlst_http()) do f
    for l in eachline(f)
      for row in split(l, "<td>")
        m = match(r"<a href='http://www.cgmlst.org/ncs/schema/(?<id>[0-9]+)/'><em>(?<species>.+)</em>", row)
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
  count = 0
  for (id, species) in available_cgmlst_schema()
    if prefix == nothing || startswith(species,prefix)
      count += 1
      @printf "%-30s - ID:%s\n" species id
    end
  end
  println("$count schema found.")
end

function _find_cgmlst_id(target_id)
  for (id,species) in available_cgmlst_schema()
    if id == target_id || species == target_id
      return id
    end
  end
  return nothing
end

function download_cgmlst_scheme(target_id, output_dir, overwrite=false)
  id = _find_cgmlst_id(target_id)
  if id == nothing
    Lumberjack.warn("Id/species ($target_id) not found!")
    exit(-1)
  end
  info("Downloading cgMLST scheme ...")
  gzip_alleles = _download_to_folder("http://www.cgmlst.org/ncs/schema/$id/alleles",output_dir)
  # unzip file to one FASTA per locus:
  loci_files = String[]
  current_fasta_fh = nothing
  locus = ""
  fh = GZip.open(gzip_alleles)
  info("Unzipping cgMLST scheme into individual FASTA files for each loci ...")
  n_locus = 0
  for l in eachline(fh)
    # skip empty lines or "="
    if length(strip(l)) == 0 || startswith(l, "=")
      continue
    end
    if startswith(l,"#")
      # print to show progress:
      if n_locus % 200 == 0
        print(".")
      end
      n_locus += 1
      locus = strip(l[2:end]) # remove the starting '#'
      if current_fasta_fh != nothing
        close(current_fasta_fh)
      end
      fasta = joinpath(output_dir, "$locus.fa")
      push!(loci_files, fasta)
      current_fasta_fh = open(fasta, "w")
    elseif startswith(l,">")
      # substitute the number only to locus_number
      id = strip(l[2:end])
      write(current_fasta_fh, ">$(locus)_$id\n")
    else #
      if current_fasta_fh != nothing
        write(current_fasta_fh, l)
      end
    end
  end
  println()
  info("$n_locus loci found.")
  return loci_files
end
