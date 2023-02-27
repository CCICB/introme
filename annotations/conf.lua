function contains(str, sub)
	str = str:gsub("%W", "")
	sub = sub:gsub("%W", "")
	return string.find(str, sub) ~= nil
end

function split(str, sep)
        local sep, fields = sep or ":", {}
        local pattern = string.format("([^%s]+)", sep)
        str:gsub(pattern, function(c) fields[#fields+1] = c end)
        return fields
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


function spliceai(entry) -- processes precomputed SpliceAI scores
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
