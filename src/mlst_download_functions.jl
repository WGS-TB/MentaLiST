using LightXML

function _get_first_line(node)
  return split(content(node), "\n")[1]
end

function _pubmlst_xml()
  if !isfile("dbases.xml")
    download("https://pubmlst.org/data/dbases.xml", "dbases.xml")
  end
  return parse_file("dbases.xml")
end

function list_pubmlst_schema(prefix)
  xroot = root(_pubmlst_xml())
  species_list = String[]
  for species in xroot["species"]
    sp_name = _get_first_line(species)
    if prefix == nothing || startswith(sp_name, prefix)
      push!(species_list,sp_name)
    end
  end
  println(join(species_list,"\n"))
  println("$(length(species_list)) schema found.")
end


function _download_to_folder(url, output_dir, overwrite=false)
  filepath = joinpath(output_dir, basename(url))
  if overwrite || !isfile(filepath)
    download(url, filepath)
  end
  return filepath
end

function download_pubmlst_scheme(target_species, output_dir, overwrite=false)
  loci_files = String[]
  xroot = root(_pubmlst_xml())
  info("Searching for the scheme ... ")
  for species in xroot["species"]
    sp_name = _get_first_line(species)
    # mlst_path = joinpath(output_dir,sp_name)
    if sp_name != target_species
      continue
    end
    info("Downloading ...")
    mkpath(output_dir)
    db = species["mlst"][1]["database"][1]
    profile_url = content(db["profiles"][1]["url"][1])
    _download_to_folder(profile_url, output_dir, overwrite)
    for locus in db["loci"][1]["locus"]
      locus_name = _get_first_line(locus)
      info("Downloading locus $locus_name ...")
      locus_url = content(locus["url"][1])
      filepath = _download_to_folder(locus_url, output_dir)
      push!(loci_files, filepath)
    end
  end
  info("Finished downloading.")
  return loci_files
end
