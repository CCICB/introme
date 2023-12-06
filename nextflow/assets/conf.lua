function contains(str, sub)
	str = str:gsub("%W", "")
	sub = sub:gsub("%W", "")
	return string.find(str, sub) ~= nil
end

function split(str, sep)
    local sep, fields = sep or ":", {}
    local pattern = string.format("([^%s]+)", sep)
	if str ~= nil and str ~= '' then
		str:gsub(pattern, function(c) fields[#fields+1] = c end)
		return fields
	end
end

function indexOf(table, value) -- returns the first index of the table that matched the value
	for i, v in ipairs(table) do
		if v == value then
			return i
		end
	end
	return nil
end

function location(str)
	local fields = split(str,",")
	for i=1,#fields do
		if contains(fields[i], "exon") then
			return "exonic"
		end
	end

	for i=1,#fields do
		if contains(fields[i], "gene") then
        	return "intronic"
        end
	end
end

function find(str, search)
	local fields = split(str, " ")
	local result = ""
	for i=1,#fields do
		if contains(fields[i], search) then
			local match = fields[i+1]:gsub("[^%w_-]","")
			if result == "" then
				result = match
			elseif not(contains(result, match)) then
				result = result .. "," .. match
			end
		end
	end
	return result
end

function regions(region_string)
	local t = type(region_string)
	if t == "string" then
		return region_string
	elseif t == "table" then
		for i=1,#region_string do
			if contains(region_string[i], "donor_canonical") then
				return "donor_canonical"
			elseif contains(region_string[i], "acceptor_canonical") then
				return "acceptor_canonical"
			end
		end
		for i=1,#region_string do
			if contains(region_string[i], "donor_exonic") then
				return "donor_exonic"
			elseif contains(region_string[i], "acceptor_exonic") then
				return "acceptor_exonic"
			end
		end
		for i=1,#region_string do
			if contains(region_string[i], "donor_region") then
				return "donor_region"
			elseif contains(region_string[i], "acceptor_region") then
				return "acceptor_region"
			end
		end
		for i=1,#region_string do
			if contains(region_string[i], "branchpoint_region") then
				return "branchpoint_region"
			end
		end
	end
end


function spliceai(entry) -- processes SpliceAI scores
	local t = type(entry)
	if t == "string" then
		return entry -- returns original value if single entry
	elseif t == "table" then
		local maximums = {}
		for i=1,#entry do -- calculate the maximum SpliceAI score of AG, AL, DG, DL for each entry
			maximums[i] = math.max(tonumber(split(entry[i], "|")[3]), tonumber(split(entry[i], "|")[4]), tonumber(split(entry[i], "|")[5]), tonumber(split(entry[i], "|")[6]))
		end
		return entry[indexOf(maximums, math.max(unpack(maximums)))] -- returns the full record which contains the maximum SpliceAI score
	end
end


function mmsplice(entry) -- processes MMSplice scores
	local t = type(entry)
	if t == "string" then
		return entry -- returns original value if single entry
	elseif t == "table" then
		local maximums = {}
		for i=1,#entry do -- calculate the entry with the maximum absolute delta_logit_PSI
			maximums[i] = math.abs(tonumber(split(entry[i], "|")[6]))
		end
		return entry[indexOf(maximums, math.max(unpack(maximums)))] -- returns the full record which contains the maximum MMSplice score
	end
end

function pangolin(entry) -- processes Pangolin scores
	local t = type(entry)
	if t == "string" then
		return entry -- returns original value if single entry
	elseif t == "table" then
		local maximums = {}
		for i=1,#entry do -- calculate the maximum Pangolin score of GAIN_POS, GAIN_SCORE, LOSS_POS, LOSS_SCORE for each entry
			local gain = split(entry[i], "|")[2]
			local loss = split(entry[i], "|")[3]
			maximums[i] = math.max(tonumber(split(gain, ":")[2]), math.abs(tonumber(split(loss, ":")[2])))
		end
		return entry[indexOf(maximums, math.max(unpack(maximums)))] -- returns the full record which contains the maximum Pangolin score
	end
end

function spip(entry) -- processes precomputed SPiP scores
	local t = type(entry)
	if t == "string" then
		return entry -- returns original value if single entry
	elseif t == "table" then
		local interconfident = tonumber(split(entry[i], "|")[4])
		maximums[i] = tonumber(split(interconfident, "%")[1]) * 1.0
		return entry[indexOf(maximums, math.max(unpack(maximums)))] -- returns the full record which contains the maximum InterConfident score
	end
end

function spip_interconfident(str) -- get the InterConfident Score
	local interconfident = tonumber(split(str, "|")[4])
	return tonumber(split(interconfident, "%")[1]) * 1.0
end

function spip_min(str) -- get the InterConfident Min Range
	local interconfident = tonumber(split(str, "|")[4])
	local segment = tonumber(split(interconfident, "%")[2])
	return tonumber(string.match(segment, "%d+%.?%d*")) * 1.0
end

function spip_max(str) -- get the InterConfident Max Range
	local interconfident = tonumber(split(str, "|")[4])
	local segment = tonumber(split(interconfident, "%")[3])
	return tonumber(string.match(segment, "%d+%.?%d*")) * 1.0
end