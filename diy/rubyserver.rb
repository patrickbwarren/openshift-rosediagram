#!/usr/bin/env ruby

require 'webrick'
include WEBrick

Dir.chdir(ARGV[1])

config = {}

config.update(:Port => 8080)
config.update(:BindAddress => ARGV[0])
config.update(:DocumentRoot => ARGV[1])

server = HTTPServer.new(config)

server.mount_proc '/exec' do |req, res|
  File.open('url.txt', 'w') do |f|
    f.puts(req.query["url"])
  end
  res['Content-Type'] = "text/plain"
  res.body = `cd #{ARGV[1]}; /usr/bin/perl wrap.pl --ns=#{ req.query["ns"] } --color=#{ req.query["color"] } --title="#{ req.query["title"] }"`
end

['INT', 'TERM'].each {|signal| trap(signal) {server.shutdown}}

server.start
