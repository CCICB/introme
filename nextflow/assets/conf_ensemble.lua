function contains(str, sub)
	str = str:gsub("%W", "")
	sub = sub:gsub("%W", "")
	return string.find(str, sub) ~= nil
end

-- Split a string by the given separator
function split(str, sep)
	local sep = sep or ":"
	local fields = {}
	local pattern = string.format("([^%s]+)", sep)
	if str ~= nil and str ~= '' then
	  str:gsub(pattern, function(c) fields[#fields+1] = c end)
	end
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


-- Given a SpliceAI annotation string, split it into transcript records.
-- This function assumes that each record has exactly `expectedFields` fields
-- (i.e. expectedFields - 1 pipe delimiters). If this assumption is violated,
-- an error will be raised.
-- The reason for this function is because sometimes there are commas for multiple
-- Ensembl transcript codes within a single spliceai entry so commas are not reliable
function splitSpliceAIRecords(str)
    local expectedFields = 19  -- Adjust this if spliceai format changes.
    local records = {}
    local pos = 1
    local len = #str

    while pos <= len do
        local start = pos
        local fieldCount = 0

        -- Scan until we have encountered (expectedFields - 1) pipes.
        while pos <= len and fieldCount < expectedFields - 1 do
            local c = string.sub(str, pos, pos)
            if c == "|" then
                fieldCount = fieldCount + 1
            end
            pos = pos + 1
        end

        if fieldCount < expectedFields - 1 then
            error(string.format("Unexpected end of string: expected at least %d pipes but found only %d starting at position %d.", expectedFields - 1, fieldCount, start))
        end

        -- Now, read until we hit a comma (which should separate transcript records)
        -- or reach the end of the string. This reads the remainder of the 19th field.
        local recordEnd = pos
        while recordEnd <= len do
            local c = string.sub(str, recordEnd, recordEnd)
            if c == "," then
                break
            end
            recordEnd = recordEnd + 1
        end

        local record = string.sub(str, start, recordEnd - 1)

        -- Check if the record splits into exactly the expected number of fields.
        local fields = split(record, "|")
        if #fields ~= expectedFields then
            error(string.format("Record does not have the expected number of fields (%d): %s", expectedFields, record))
        end

        table.insert(records, record)
        pos = recordEnd + 1  -- Move past the comma (if present)
		-- Debug print statement here:
		print(string.format("Parsed record: %s", record))
    end

    return records
end

-- Process the SpliceAI field
function spliceai(entry)
    -- If vcfanno pre-split the annotation into a table, join it back together.
    if type(entry) == "table" then
        entry = table.concat(entry, ",")
    end

    -- Now entry is a string. If it contains a comma, it may have multiple records.
    if type(entry) == "string" then
        if string.find(entry, ",") then
            entry = splitSpliceAIRecords(entry)
        else
            return entry
        end
    end

  -- Now entry should be a table of transcript strings
	if type(entry) == "table" then
		local best = nil
    	local best_score = -math.huge  -- start with a very low score
    	for i=1, #entry do
			local fields = split(entry[i], "|")
			local ds_ag = tonumber(fields[4])
			local ds_al = tonumber(fields[5])
			local ds_dg = tonumber(fields[6])
			local ds_dl = tonumber(fields[7])

			if ds_ag == nil or ds_al == nil or ds_dg == nil or ds_dl == nil then
				error("One or more DS values are missing in annotation: " .. entry[i])
			end

			local maximum = math.max(ds_ag, ds_al, ds_dg, ds_dl)
			if maximum > best_score then
				best_score = maximum
				best = i
			end
    	end
		-- Find the annotation with the highest DS score.
		return entry[best]
	else
		error("was expecting entry to be a table: " .. entry)
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
