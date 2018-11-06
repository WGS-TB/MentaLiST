# Packages needed for parallel processing need to be loaded "twice":s
using Bio.Seq: DNASequence, DNAKmer, canonical, FASTASeqRecord, FASTAReader, each, neighbors
@everywhere using Bio.Seq: DNASequence, DNAKmer, canonical, FASTASeqRecord, FASTAReader, each, neighbors
using DataStructures: DefaultDict
@everywhere using DataStructures: DefaultDict
using FastaIO
@everywhere using FastaIO

@everywhere function twin{k}(::Type{DNAKmer{k}}, km)
    DNAKmer{k}(Bio.Seq.reverse_complement(DNASequence("$km")))
end

function contig_to_string(c)
  return "$(c[1])" * join(["$(x[end])" for x in c[2:end]])
end

@everywhere function get_contig{k}(::Type{DNAKmer{k}}, kmers, km)
  contig_fw = get_contig_forward(DNAKmer{k}, kmers, km)
  contig_bw = get_contig_forward(DNAKmer{k}, kmers, twin(DNAKmer{k}, km))

  if km in neighbors(contig_fw[end])
    return contig_fw
  else
    return [[twin(DNAKmer{k}, x) for x in reverse(contig_bw)[1:end-1]] ; contig_fw]
  end
end

@everywhere function get_contig_forward{k}(::Type{DNAKmer{k}}, kmers, km)
    c_fw = DNAKmer{k}[km]
    while true
        # check if forward exists, and only 1:
        fw_neighbors = [kmer for kmer in neighbors(c_fw[end]) if canonical(kmer) in kmers]
        if length(fw_neighbors) != 1
          break
        end
        candidate = fw_neighbors[1]
        if candidate == km || candidate == twin(DNAKmer{k}, km) || candidate == twin(DNAKmer{k}, c_fw[end])
          break
        end
        bw_neighbors = [kmer for kmer in neighbors(twin(DNAKmer{k}, candidate)) if canonical(kmer) in kmers]
        if length(bw_neighbors) != 1
          break
        end
        push!(c_fw, candidate)
    end
    return c_fw
end

@everywhere function all_contigs{k}(::Type{DNAKmer{k}}, kmers)
  done = Set{DNAKmer{k}}()
  contigs = []
  for kmer in kmers
    if kmer in done
      continue
    end
    contig = get_contig(DNAKmer{k}, kmers, kmer)
      for c_kmer in contig
      push!(done, canonical(c_kmer))
    end
    push!(contigs, contig)
  end
  return contigs
end

### classifies each kmer with his colors (present alleles)
@everywhere function find_all_kmer_colors{k}(::Type{DNAKmer{k}}, fastafile)
  record = FASTASeqRecord{DNASequence}()
  kmer_class = DefaultDict{DNAKmer{k}, Vector{Int16}}(() -> Int16[])
  allele_ids = Int16[]
  allele_idx::Int16 = 1
  open(FASTAReader{DNASequence}, fastafile) do reader
    while !eof(reader)
      try
        read!(reader, record)
        seen = Set{DNAKmer{k}}()
        for (pos, kmer) in each(DNAKmer{k}, record.seq)
          can_kmer = canonical(kmer)
          if !in(can_kmer, seen)
            push!(kmer_class[can_kmer], allele_idx)
            push!(seen, can_kmer)
          end
        end
      catch
        println("CRITICAL ERROR: Error parsing file $fastafile, at record $(record.name), most likely some unkown characters. Please fix it and try again.")
        exit(-1)
      end
      # find the separator; will assume that if I see a "_", that's it, otherwise try "-";
      separator = in('_', record.name) ? "_" : "-"
      # update idx; the counter idx is incremental (1,2, ...) because we need the array sorted.
      # But this is not always sin the allele ordering, so we have to save the original id to restore it later;
      allele_id = parse(Int16,split(record.name, separator)[end])
      push!(allele_ids, allele_id)
      allele_idx += 1
    end
  end
  return kmer_class, allele_ids
end

@everywhere function build_db_graph{k}(::Type{DNAKmer{k}}, fastafile)
  # get kmers and colors (alleles)
  kmer_class, allele_ids = find_all_kmer_colors(DNAKmer{k}, fastafile)

  contig_list = all_contigs(DNAKmer{k}, keys(kmer_class))
  # select 1 kmer per contig, and set weights:
  kmer_weights = Dict{UInt64, Int}() # Kmer is converted to UInt64, because DNAKmer apparently does not support write(), needed for the parallel pmap() call;
  filtered_kmer_class = Dict{UInt64, Vector{Int16}}()
  for contig in contig_list
      kmer = canonical(contig[1])
      kmer_uint = convert(UInt64,kmer)
      kmer_weights[kmer_uint] = length(contig)
      filtered_kmer_class[kmer_uint] = kmer_class[kmer]
  end
  return (filtered_kmer_class, allele_ids, kmer_weights)
end
