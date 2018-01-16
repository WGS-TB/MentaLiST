function twin{k}(::Type{DNAKmer{k}}, km)
    DNAKmer{k}(Bio.Seq.reverse_complement(DNASequence("$km")))
end

function build{k}(::Type{DNAKmer{k}}, fastafiles)
  kmers = Set{DNAKmer{k}}()
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  for fastafile in fastafiles
    open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
      while !eof(reader)
        read!(reader, record)
        for (pos, kmer) in each(DNAKmer{k}, record.seq)
          push!(kmers, kmer)
        end
        for (pos, kmer) in each(DNAKmer{k}, Bio.Seq.reverse_complement(record.seq))
          push!(kmers, kmer)
        end
      end
    end
  end
  return kmers
end

function contig_to_string(c)
  return "$(c[1])" * join(["$(x[end])" for x in c[2:end]])
end

function get_contig{k}(::Type{DNAKmer{k}}, kmers, km)
  contig_fw = get_contig_forward(DNAKmer{k}, kmers, km)
  contig_bw = get_contig_forward(DNAKmer{k}, kmers, twin(DNAKmer{k}, km))

  if km in neighbors(contig_fw[end])
    return contig_fw
  else
    return [[twin(DNAKmer{k}, x) for x in reverse(contig_bw)[1:end-1]] ; contig_fw]
  end
end

function get_contig_forward{k}(::Type{DNAKmer{k}}, kmers, km)
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

function all_contigs{k}(::Type{DNAKmer{k}}, kmers)
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
function find_all_kmer_colors{k}(::Type{DNAKmer{k}}, fastafile)
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  kmer_class = DefaultDict{DNAKmer{k}, Vector{Int16}}(() -> Int16[])
  allele_ids = Int16[]
  allele_idx::Int16 = 1
  open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
    while !eof(reader)
      read!(reader, record)
      seen = Set{DNAKmer{k}}()
      for (pos, kmer) in each(DNAKmer{k}, record.seq)
        can_kmer = canonical(kmer)
        if !in(can_kmer, seen)
          push!(kmer_class[can_kmer], allele_idx)
          push!(seen, can_kmer)
        end
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


function build_db_graph{k}(::Type{DNAKmer{k}}, fastafile)
  # get kmers and colors (alleles)
  kmer_class, allele_ids = find_all_kmer_colors(DNAKmer{k}, fastafile)

  contig_list = all_contigs(DNAKmer{k}, keys(kmer_class))
  # select kmers:
  kmer_weights = Dict{DNAKmer{k}, Int}()
  # filter kmer class to only kmers that are selected per contig;
  filtered_kmer_class = Dict{DNAKmer{k}, Vector{Int16}}()
  for contig in contig_list
    # select 1st kmer of contig:
    can_kmer = canonical(contig[1])
    kmer_weights[can_kmer] = length(contig)
    filtered_kmer_class[can_kmer] = kmer_class[can_kmer]
  end
  return filtered_kmer_class, allele_ids, kmer_weights
end

# Calling:
function sequence_coverage{k}(::Type{DNAKmer{k}}, sequence, kmer_count, kmer_thr=6)
  # TODO: not using the gap list anymore, just for debug. Eventually delete!

  # find the minimum coverage depth of all kmers of the given variant, and
  # how many kmers are covered/uncovered
  gap_list = Tuple{Int,Int}[]
  smallest_coverage = 100000
  covered = 0
  uncovered = 0

  in_gap = false
  gap_start = 0

  for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
    can_kmer = canonical(kmer)
    cov = get(kmer_count, can_kmer, 0)
    if cov >= kmer_thr
      covered += 1
      # gap:
      if in_gap
        push!(gap_list, (gap_start, pos-1))
        in_gap = false
      end
    else
      uncovered += 1
      # gap
      if !in_gap
        gap_start = pos
        in_gap = true
      end
    end
    if cov < smallest_coverage
      smallest_coverage = cov
    end
  end
  if in_gap
    push!(gap_list, (gap_start, length(sequence)-k+1))
  end
  # merge too close gaps: if by chance we get a kmer match in the middle of a gap, this will break a gap in two, even
  # with just one mutation, giving errors. TODO: How to prevent that?
  # simple idea here, merge gaps (s1,e1) and (s2,e2) if s1 + k >= s2;
  merged_gap_list = []
  i = 1
  while i <= length(gap_list)
    # if too close to the next, merge:
    if i<length(gap_list) && gap_list[i][1] + k >= gap_list[i+1][1]
      push!(merged_gap_list, (gap_list[i][1], gap_list[i+1][2]))
      i += 2
    else
      push!(merged_gap_list, gap_list[i])
      i += 1
    end
  end
  return smallest_coverage, covered, uncovered, merged_gap_list
end

function cover_sequence_gap{k}(::Type{DNAKmer{k}}, sequence, kmer_count, kmer_thr=8, max_mutations=4)
  gap_cover_list = []
  found = Set{String}()
  function scan_variant(n_mut, sequence, events, start_pos=1)
    if n_mut > max_mutations
      return
    end
    last_kmer_present = false
    for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
      if pos < start_pos
        continue
      end
      can_kmer = canonical(kmer)
      if get(kmer_count, can_kmer, 0) >= kmer_thr
        if pos > start_pos +1 && (!last_kmer_present)
          # found a kmer and the last was not there; check if I can find kmer variations for the previous position
          km = "$kmer"
          for base in ["A","C","T","G"]
            if get(kmer_count, canonical(DNAKmer{k}("$base$(km[1:end-1])")), 0) >= kmer_thr
              # found a kmer variant;
              # Test one more kmer before adding to the pool:
              # substitution at position pos-1
              substitution_seq = sequence[1:pos-2] * base * sequence[pos:end]
              push!(candidate_list, (n_mut+1, substitution_seq, [events; [("S",pos-1,"$(sequence[pos-1])->$base")]], start_pos))
              # insertion:
              insertion_seq = sequence[1:pos-1] * base * sequence[pos:end]
              push!(candidate_list, (n_mut+1, insertion_seq, [events; [("I",pos-1,base)]], start_pos))
              # deletion:
              for i in 0:2
                if pos-2-i < 1
                  break
                end
                if sequence[pos-2-i:pos-2-i] == base
                  deletion_seq = sequence[1:pos-2-i] * sequence[pos:end]
                  push!(candidate_list, (n_mut+1+i, deletion_seq, [events; [("D",pos-1,"$(i+1)")]], start_pos))
                  break
                end
              end
              # this variant has mismatches, so just return;
              return
            end
          end
        end
        last_kmer_present = true
      else
        if last_kmer_present
          # last_kmer_present was present, but this not; check kmer variations for nest position
          km = "$kmer"
          for base in ["A","C","T","G"]
            if get(kmer_count, canonical(DNAKmer{k}("$(km[1:end-1])$base")), 0) >= kmer_thr
              # println("Found a $base at pos $(pos+k-1)")
              # substitution at position pos+k-1;
              substitution_seq = sequence[1:pos+k-2] * base * sequence[pos+k:end]
              push!(candidate_list, (n_mut+1, substitution_seq, [events; [("S",pos+k-1,"$(sequence[pos+k-1])->$base")]], pos-1))
              # insertion:
              insertion_seq = sequence[1:pos+k-2] * base * sequence[pos+k-1:end]
              push!(candidate_list, (n_mut+1, insertion_seq, [events; [("I",pos+k-1,base)]], pos-1))
              # deletion:
              for i in 0:2
                # println("trying $(pos+k+i), found $(sequence[pos+k+i]) (want $base) $(sequence[pos+k+i:pos+k+i] == base).")
                if pos+k+i > length(sequence)
                  break
                end
                if sequence[pos+k+i:pos+k+i] == base
                  # println("Found a possible deletion, pos $(pos+k-1) to $(pos+k+i-1)")
                  deletion_seq = sequence[1:pos+k-2] * sequence[pos+k+i:end]
                  push!(candidate_list, (n_mut+1+i, deletion_seq, [events; [("D",pos+k-1,"$(i+1)")]], pos-1))
                  break
                end
              end
              # the original sequence has mismatches, so just return;
              return
            end
          end
        end
        last_kmer_present = false
      end
    end
    # passed without suggesting mutations, test this variant:
    seq_cov = sequence_coverage(DNAKmer{k}, sequence, kmer_count)[1]
    if seq_cov[1] >= kmer_thr
      if !in(sequence, found)
        push!(gap_cover_list, (n_mut, sequence, events, seq_cov))
        push!(found, sequence)
        max_mutations = n_mut
      end
    end
  end
  candidate_list = [(0, sequence, [], 0)]
  idx = 1
  while idx <= length(candidate_list)
    n_mut, candidate_sequence, mutations, start_pos = candidate_list[idx]
    scan_variant(n_mut, candidate_sequence, mutations, start_pos)
    idx += 1
  end

  if length(gap_cover_list) > 0
    # sort and return the one with minimum mutations: TODO: tie?
    return sort(gap_cover_list, by=x->x[1])[1]
  else
    return nothing
  end

end
# struct CorrectedTemplate
#
# end
function correct_template{k}(::Type{DNAKmer{k}}, template_seq, gap_list, kmer_count, kmer_thr, max_mutations)
  function find_next_gap(sequence, skip=1)
    in_gap = false
    gap_start = 0
    for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
      if pos < skip
        continue
      end
      can_kmer = canonical(kmer)
      cov = get(kmer_count, can_kmer, 0)
      if cov >= kmer_thr
        if in_gap
          return (gap_start, pos-1)
        end
      else
        if !in_gap
          gap_start = pos
          in_gap = true
        end
      end
    end # loop
    if in_gap
      return (gap_start, length(sequence)-k+1)
    end
    return nothing
  end # end auxiliar function find_next_gap

  # start correction:
  current_skip = 1
  uncorrected_gaps = []
  corrected_seq = template_seq
  total_mut = 0
  mutations_list = []
  while current_skip < length(corrected_seq)
    gap = find_next_gap(corrected_seq, current_skip)
    if gap == nothing
      break
    end
    st_pos, end_pos = gap
    adj_start = max(st_pos-1, 1)
    adj_end = min(end_pos+k, length(corrected_seq))
    gap_seq = corrected_seq[adj_start:adj_end]
    gap_cover = cover_sequence_gap(DNAKmer{k}, gap_seq, kmer_count, kmer_thr, max_mutations)
    if gap_cover == nothing
      push!(uncorrected_gaps, (st_pos, end_pos))
      current_skip = end_pos + 1
    else
      # TODO: report n_mut, mut_list, depth, etc.
      n_mut, gap_cover_seq, mut_list, depth = gap_cover
      # println("Mutations found: $mut_list")
      total_mut += n_mut
      for mut in mut_list
        # TODO: correct idx
        push!(mutations_list, (mut[1], mut[2] + adj_start - 1, mut[3]))
      end
      corrected_seq = corrected_seq[1:adj_start-1] * gap_cover_seq * corrected_seq[adj_end+1:end]
      current_skip = adj_start + length(gap_cover_seq) - k
    end
  end
  return corrected_seq, total_mut, mutations_list, uncorrected_gaps
end

function describe_mutation(mut)
  mut, pos, desc = mut
  if mut == "D"
    return "Del of len $desc at pos $pos"
  elseif mut == "I"
    return "Ins of base $desc at pos $pos"
  elseif mut == "S"
    return "Subst $desc at pos $pos"
  else
    return ""
  end
end
