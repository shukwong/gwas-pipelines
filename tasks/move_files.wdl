version development

workflow move_outputs {
    input {
        File rclone_box_config
        File input_file
        String box_directory 
    }

    call move_file_to_box {
      input:
        rclone_box_config = rclone_box_config,
        input_file = input_file,
        box_directory = box_directory
    }

     output {
        String file_location_on_box = move_file_to_box.file_location_on_box
    }
}    

task move_file_to_box {
    
    input {
      File rclone_box_config
      File input_file
      String box_directory
    }
    

    command <<<
        rclone --config=~{rclone_box_config} copy ~{input_file} box:~{box_directory}
        rclone --config=~{rclone_box_config} ls box:~{box_directory} > file_location_on_box.txt
    >>>    

	runtime {
		docker: "quay.io/shukwong/rclone"
	}

  output {
	    String file_location_on_box = "file_location_on_box.txt"
  }
}

